#!/usr/bin/env python3
"""Takes in any set of timestamp data, and reads/converts it.

Implemented for timestamp7 data, i.e. 4ps resolution data, but also works for other
timestamps, e.g. 2ns and 125ps resolutions. Note the timestamp formats:
'-a0' hex text, '-a1' binary, '-a2' binary text.

Usage:

  Split timestamp data on same channel into two different timestamp files:

    >>> t, p = read_a1("timestampfile.Aa1X.dat", legacy=True)

    >>> mask = (p & 1).astype(bool)     # events with bit 0 set
    >>> write_a1("alice.Aa1.dat", t[mask], p[mask])

    >>> mask = p == 8                   # events with *only* bit 3 set
    >>> mask &= ((t > 100) & (t < 120)  # events with timestamp between 100s and 120s
    >>> write_a2("bob.Aa2X.txt", t[mask], p[mask], legacy=True)


  Equivalent command-line usage (faster since using smaller stream buffers):

    > ./parse_timestamps.py -A1 -X -a1 --pattern 1 --mask \
         timestampfile.Aa1X.dat alice.Aa1.dat

    > ./parse_timestamps.py -A1 -X -a1 -x --pattern 8 --start 100 --end 120 \
         timestampfile.Aa1X.dat bob.Aa2X.txt


  Read large timestamp files in Python (streams in small buffers to avoid OOM):

    >>> stream, num_batches = sread_a1("timestampfile.Aa1X.dat", legacy=True)
    >>> for i, (t, p) in enumerate(stream, start=1):
    ...     print(f"Batch: {i}/{num_batches}")
    ...     # Do stuff with t and p


  Print event statistics:

    > ./parse_timestamps.py -X timestampfile.Aa1X.dat


Changelog:
    2022-10-07 Justin: Dummy events and blinding events all set to 0, not handled
    2023-05-09 Justin: Add streaming capabilities, minimize inline processing
    2024-07-09 Justin: Remove references to deprecated typing.Tuple types
"""

import argparse
import bisect
import io
import pathlib
import struct
import sys
import warnings
from typing import Optional

import numba
import numpy as np
import tqdm

from fpfind import NP_PRECISEFLOAT, TSRES

###########################################
#  PARSERS (from internal integer array)  #
###########################################


def _parse_a1(
    data,
    legacy: bool = False,
    resolution: TSRES = TSRES.NS1,
    fractional: bool = True,
    ignore_rollover: bool = False,
):
    """
    Args:
        data: np.array[np.int32] with shape (N, 2)
    """
    high_pos = 1
    low_pos = 0
    if legacy:
        high_pos, low_pos = low_pos, high_pos

    if ignore_rollover:
        r = (data[:, low_pos] & 0b10000).astype(bool)  # check rollover flag
        data = data[~r]

    t = (np.uint64(data[:, high_pos]) << 22) + (data[:, low_pos] >> 10)
    t = _format_timestamps(t, resolution, fractional)
    p = data[:, low_pos] & 0xF
    return t, p


def _parse_a0(
    data,
    resolution: TSRES = TSRES.NS1,
    fractional: bool = True,
    ignore_rollover: bool = False,
):
    """
    Args:
        data: np.array[np.int32] with shape (N, 2)
    """
    return _parse_a1(data, False, resolution, fractional, ignore_rollover)


def _parse_a2(
    data,
    resolution: TSRES = TSRES.NS1,
    fractional: bool = True,
    ignore_rollover: bool = False,
):
    """
    Args:
        data: np.array[np.int64] with shape (N,)
    """
    if ignore_rollover:
        r = (data & 0b10000).astype(bool)  # check rollover flag
        data = data[~r]

    t = np.uint64(data >> 10)
    t = _format_timestamps(t, resolution, fractional)
    p = data & 0xF
    return t, p


def _format_timestamps(t: list, resolution: TSRES, fractional: bool):
    """Returns conversion of timestamps into desired format and resolution.

    Note:
        Short-circuit when raw timestamps are requested is applied, i.e.
        'resolution' of TSRES.PS4 and 'fractional' is False. Avoiding computation
        across the entire array has significant time-savings, see below:

        >>> import timeit
        >>> t = np.random.randint(0, int(1e18), size=(100_000,), dtype=np.uint64)
        >>> f = lambda: _format_timestamps(t, TSRES.PS4, False)
        >>> _ = timeit.timeit(f, number=10_000)
        # 430 us per loop without short-circuit
        # 989 ns per loop with short-circuit

        In practice, the number of batches to process is not particularly high,
        e.g. 1GB timestamp file corresponds to 1250 loops.
    """
    if fractional:
        t = np.array(t, dtype=NP_PRECISEFLOAT)
        t = t / (TSRES.PS4.value / resolution.value)
    elif resolution is not TSRES.PS4:  # short-circuit
        t = np.array(t, dtype=np.uint64)
        t = t // int(TSRES.PS4.value // resolution.value)
    return t


########################
#  UNBUFFERED READERS  #
########################


def read_a0(
    filename: str,
    legacy: bool = None,
    resolution: TSRES = TSRES.NS1,
    fractional: bool = True,
    ignore_rollover: bool = False,
):
    """Converts a0 timestamp format into timestamps and detector pattern.

    'legacy' is set at the point of 'readevents' invocation, and determines
    the timestamp storage format.

    'resolution' can be set to any of the available TSRES enumeration, and
    denotes the units of the output timestamps. 'fractional' determines
    whether sub-unit precisions should be preserved as well. The combination
    of 'resolution' and 'fractional' will determine the resulting output
    format of the timestamps.

    Note:
        128-bit floating point is required since 64-bit float only has a
        precision of 53-bits (timestamps have precision of 54-bits). If
        the bandwidth cost is too high when using 128-bit floats, then
        chances are the application also does not require sub-ns precision.

        There is no legacy formatting for event types a0 and a2. Preserved in
        the function signature to maintain parity with 'read_a1'.
    """
    data = np.genfromtxt(filename, delimiter="\n", dtype="U8")
    data = np.array([int(v, 16) for v in data]).reshape(-1, 2)
    return _parse_a0(data, resolution, fractional, ignore_rollover)


def _read_a1(
    filename: str,
    legacy: bool = False,
    resolution: TSRES = TSRES.NS1,
    fractional: bool = True,
    ignore_rollover: bool = False,
):
    """See documentation for 'read_a0'."""
    with open(filename, "rb") as f:
        data = np.fromfile(file=f, dtype="=I").reshape(-1, 2)
    return _parse_a1(data, legacy, resolution, fractional, ignore_rollover)


def read_a2(
    filename: str,
    legacy: bool = None,
    resolution: TSRES = TSRES.NS1,
    fractional: bool = True,
    ignore_rollover: bool = False,
):
    """See documentation for 'read_a0'."""
    data = np.genfromtxt(filename, delimiter="\n", dtype="U16")
    data = np.array([int(v, 16) for v in data])
    return _parse_a2(data, resolution, fractional, ignore_rollover)


######################
#  BUFFERED READERS  #
######################


def read_a1(
    filename: str,
    legacy: bool = False,
    resolution: TSRES = TSRES.NS1,
    fractional: bool = True,
    start: float = None,
    end: float = None,
    pattern: int = None,
    mask: bool = False,
    invert: bool = False,
    buffer_size: int = 100_000,
    relative_time: bool = False,
    ignore_rollover: bool = False,
):
    """Wrapper to sread_a1 for efficiently reading subsets of whole file.

    For descriptions of arguments, read the documentation:

        - 'read_a0': Importing timestamps.
        - 'sread_a1': Importing timestamps via file streaming.
        - 'get_timing_bounds': Filtering by timing.
        - 'get_pattern_mask': Filtering by detector pattern.

    If 'buffer_size' of None is supplied, the old legacy '_read_a1' code
    will be used, and all timing/pattern filtering will be disabled.

    If 'relative_time' is True, 'start' and 'end' values (in seconds) will
    be shifted relative to the first timestamp.
    """
    if buffer_size is None:
        return _read_a1(filename, legacy, resolution, fractional, ignore_rollover)

    streamer, _ = sread_a1(
        filename, legacy, resolution, fractional, buffer_size, ignore_rollover
    )
    ts = []
    ps = []
    commenced_reading = False
    has_applied_time_offset = False
    for t, p in streamer:
        if relative_time and not has_applied_time_offset:
            has_applied_time_offset = True
            if start is not None:
                start += t[0] * 1e-9
            if end is not None:
                end += t[0] * 1e-9
        i, j = get_timing_bounds(t, start, end, resolution)
        t = t[i:j]
        p = p[i:j]
        if pattern is not None:
            pmask = get_pattern_mask(p, pattern, mask, invert)
            t = t[pmask]
            p = p[pmask]

        # Terminate if no more timings readable
        if commenced_reading and len(t) == 0:
            break
        elif len(t) != 0:
            commenced_reading = True
        ts.append(t)
        ps.append(p)

    # Return array of correct dtype if no data
    if len(ts) == 0:
        t = _format_timestamps([], resolution, fractional)
        p = np.array([], dtype=np.int32)  # consistent with %I format code
        return t, p

    t = np.concatenate(ts)
    p = np.concatenate(ps)
    return t, p


def read_a1_from_buffer(
    buffer,
    legacy: bool = False,
    resolution: TSRES = TSRES.NS1,
    fractional: bool = True,
    ignore_rollover: bool = False,
):
    data = np.frombuffer(buffer, dtype="=I").reshape(-1, 2)
    return _parse_a1(data, legacy, resolution, fractional, ignore_rollover)


def read_a0_from_buffer(
    buffer,
    legacy: bool = None,
    resolution: TSRES = TSRES.NS1,
    fractional: bool = True,
    ignore_rollover: bool = False,
):
    data = buffer.strip().split("\n")
    data = np.array([int(v, 16) for v in data]).reshape(-1, 2)
    return _parse_a0(data, resolution, fractional, ignore_rollover)


def read_a2_from_buffer(
    buffer,
    legacy: bool = None,
    resolution: TSRES = TSRES.NS1,
    fractional: bool = True,
    ignore_rollover: bool = False,
):
    data = buffer.strip().split("\n")
    data = np.array([int(v, 16) for v in data])
    return _parse_a2(data, resolution, fractional, ignore_rollover)


########################
#  BUFFERED STREAMERS  #
########################


def sread_a1(
    filename: str,
    legacy: bool = False,
    resolution: TSRES = TSRES.NS1,
    fractional: bool = True,
    buffer_size: int = 100_000,
    ignore_rollover: bool = False,
):
    """Block streaming variant of 'read_a1'.

    For large timestamp datasets where either not all timestamps need
    to be loaded into memory, or only statistics need to be retrieved.
    'buffer_size' in bytes, other arguments to follow as documented in
    'read_a0'.

    If timestamp formatting is not required, the raw timestamp resolution
    can be used to speed up computation time, i.e. 'resolution' of
    TSRES.PS4 and 'fractional' of False.

    Note that these streamer can easily be extended to perform some
    pre-processing cleanup, see Usage.

    Args:
        buffer_size: Buffer size in number of events.

    Examples:
        >>> for t, p in sread_a1(...):
        ...     print(t, p)

    Note:
        Performance of naive streaming is poor (roughly 2 orders magnitude
        slower), and deprecated thusly. A compromise, whereby blocks are
        streamed instead, yields minimal performance loss while avoiding
        issues with memory (to avoid an OOM kill).

        A buffer size between 100kB and 1MB is optimal. Disk read latencies
        dominate at low buffer sizes, while memory allocation latencies
        dominate at high buffer sizes.

        Where efficiency is desired and number of timestamps is small,
        'read_a1' should be preferred instead.
    """

    def _sread_a1():
        with open(filename, "rb") as f:
            while True:
                buffer = f.read(buffer_size * 8)  # 8 bytes per event
                if len(buffer) == 0:
                    break
                yield read_a1_from_buffer(
                    buffer, legacy, resolution, fractional, ignore_rollover
                )

    size = pathlib.Path(filename).stat().st_size
    num_batches = int(((size - 1) // (buffer_size * 8)) + 1)
    return _sread_a1(), num_batches


def sread_a0(
    filename: str,
    legacy: bool = None,
    resolution: TSRES = TSRES.NS1,
    fractional: bool = True,
    buffer_size: int = 100_000,
    ignore_rollover: bool = False,
):
    """See documentation for 'sread_a1'"""

    def _sread_a0():
        with open(filename, "r") as f:
            while True:
                buffer = f.read(buffer_size * 18)  # 16 char per event + 2 newlines
                if len(buffer) == 0:
                    break
                yield read_a0_from_buffer(
                    buffer, legacy, resolution, fractional, ignore_rollover
                )

    size = pathlib.Path(filename).stat().st_size
    num_batches = int(((size - 1) // (buffer_size * 18)) + 1)
    return _sread_a0(), num_batches


def sread_a2(
    filename: str,
    legacy: bool = None,
    resolution: TSRES = TSRES.NS1,
    fractional: bool = True,
    buffer_size: int = 100_000,
    ignore_rollover: bool = False,
):
    """See documentation for 'sread_a1'"""

    def _sread_a2():
        with open(filename, "r") as f:
            while True:
                buffer = f.read(buffer_size * 17)  # 16 char per event + 1 char newline
                if len(buffer) == 0:
                    break
                yield read_a2_from_buffer(
                    buffer, legacy, resolution, fractional, ignore_rollover
                )

    size = pathlib.Path(filename).stat().st_size
    num_batches = int(((size - 1) // (buffer_size * 17)) + 1)
    return _sread_a2(), num_batches


#########################
#  AUXILIARY FUNCTIONS  #
#########################


def read_a1_start_end(
    filename: str,
    legacy: bool = False,
    resolution: TSRES = TSRES.NS1,
    fractional: bool = True,
    ignore_rollover: bool = False,
):
    t, _ = read_a1_kth_timestamp(
        filename,
        [0, -1],
        legacy,
        resolution,
        fractional,
        ignore_rollover,
    )
    return t * 1e-9


def read_a1_kth_timestamp(
    filename,
    k,
    legacy: bool = False,
    resolution: TSRES = TSRES.NS1,
    fractional: bool = True,
    ignore_rollover: bool = False,
):
    filesize = pathlib.Path(filename).stat().st_size
    num_events = filesize / 8  # 8 bytes per event

    # Parse indices
    ks = np.array(k).astype(int)
    ks = np.array(
        [(num_events + k if k < 0 else k) for k in ks]
    )  # convert to positive ints
    if np.any(ks >= num_events):
        raise ValueError(f"Index out-of-bounds, only {num_events} events available.")

    buffer = bytearray()
    with open(filename, "rb") as f:
        for _k in ks:
            f.seek(int(_k) * 8, io.SEEK_SET)  # should be O(1) in modern filesystems
            buffer.extend(f.read(8))

    return read_a1_from_buffer(buffer, legacy, resolution, fractional, ignore_rollover)


def read_a1_overlapping(
    *filenames,
    legacy: bool = False,
    resolution: TSRES = TSRES.NS1,
    fractional: bool = True,
    duration: float = None,
):
    """Reads multiple timestamp files with overlapping durations.

    Convenience function that replaces a common workflow. Streaming is used,
    and no filtering patterns are used. No detector patterns are also returned.

    If 'duration' is unspecified, the full overlapping timestamps are extracted.
    Otherwise, the timestamps will be truncated by the corresponding duration,
    in units of seconds.
    """
    # Search for start and end timings
    timings = np.array([read_a1_start_end(f, legacy=legacy) for f in filenames])
    start = np.max(timings[:, 0])
    end = np.min(timings[:, 1])
    if duration is not None:
        if duration > end - start:
            warnings.warn(
                "Duration requested is longer than available data: "
                f"{duration:.2f}s > {end-start:.2f}s"
            )
        end = duration + start

    # Extract timestamps only
    data = [
        read_a1(
            filename,
            legacy=legacy,
            resolution=resolution,
            fractional=fractional,
            start=start,
            end=end,
        )
        for filename in filenames
    ]
    timestamps, channels = list(zip(*data))
    return timestamps, channels


########################
#  UNBUFFERED WRITERS  #
########################


@numba.njit(parallel=True)
def _duplicate_event_removal(data):
    """Subroutine for _consolidate_events().

    On top of removing duplicate timestamps, this function also merges duplicate
    timestamps across different channels.

    Note:
        The numba optimization is critical: the runtimes are 31.7s (without numba),
        0.5s (with numba 1s compile), 0.4s (with numba parallel 1.6s compile), for
        a 15 Mevents timestamp file. This function is compiled lazily, i.e. only
        when duplicate removal is required.
    """
    # Identify unique timestamps and their indices
    ts = data >> np.uint64(10)
    idxs = np.flatnonzero(ts[1:] != ts[:-1]) + 1
    idxs = np.concatenate((np.array([np.int64(0)]), idxs, np.array([ts.size])))
    _t = ts[idxs[:-1]]
    _p = np.zeros_like(_t)

    # 'p' in principle can be supplied directly from '_consolidate_events()', but
    # we re-extract here from 'data' for clarity + usability. Negligible perf impact.
    p = data & np.uint64(0b1111)

    # Combine all detector patterns for the same timestamp
    for i in numba.prange(len(idxs) - 1):
        left, right = idxs[i], idxs[i + 1]
        _pi = p[left]
        for j in range(left + 1, right):
            _pi |= p[j]
        _p[i] = _pi

    # Recombine
    data = (_t << np.uint64(10)) + _p
    return data


def _consolidate_events(
    t: list,
    p: list,
    resolution: TSRES = TSRES.NS1,
    sort: bool = False,
    allow_duplicates: bool = True,
):
    """Packs events into standard a1 timestamp format.

    If sorting is required, such as during timestamp merging or timestamp
    perturbation, set 'sort' to True. This invokes O(NlogN) quick sort.

    Args:
        t: Timestamp array, in units of TSRES.NS1.
        p: Detector pattern array.
        resolution: Resolution of timestamps in timestamp array.
        sort: If events should be further sorted in chronological order.
        allow_duplicates: If events with the same timestamp should be combined,
            i.e. similar to the output from an actual time tagger.

    Note:
        float128 is needed, since float64 only encodes 53-bits of precision,
        while the high resolution timestamp has 54-bits precision.
    """
    data = (
        np.array(t, dtype=NP_PRECISEFLOAT) * (TSRES.PS4.value // resolution.value)
    ).astype(np.uint64) << 10
    data += np.array(p).astype(np.uint64)
    if sort:
        data = np.sort(data)

    # Duplicate removal can only be executed on a sorted array
    if not allow_duplicates:
        data = _duplicate_event_removal(data)

    return data


def write_a2(
    filename: str,
    t: list,
    p: list,
    legacy: bool = None,
    resolution: TSRES = TSRES.NS1,
    allow_duplicates: bool = True,
):
    data = _consolidate_events(t, p, resolution, allow_duplicates=allow_duplicates)
    with open(filename, "w") as f:
        for line in data:
            f.write(f"{line:016x}\n")


def write_a0(
    filename: str,
    t: list,
    p: list,
    legacy: bool = None,
    resolution: TSRES = TSRES.NS1,
    allow_duplicates: bool = True,
):
    events = _consolidate_events(t, p, resolution, allow_duplicates=allow_duplicates)
    data = np.empty((2 * events.size,), dtype=np.uint32)
    data[0::2] = events & 0xFFFFFFFF
    data[1::2] = events >> 32
    with open(filename, "w") as f:
        for line in data:
            f.write(f"{line:08x}\n")


def write_a1(
    filename: str,
    t: list,
    p: list,
    legacy: bool = False,
    resolution: TSRES = TSRES.NS1,
    allow_duplicates: bool = True,
):
    events = _consolidate_events(t, p, resolution, allow_duplicates=allow_duplicates)
    with open(filename, "wb") as f:
        for line in events:
            if legacy:
                line = int(line)
                line = ((line & 0xFFFFFFFF) << 32) + (line >> 32)
            f.write(struct.pack("=Q", line))


######################
#  BUFFERED WRITERS  #
######################


def swrite_a1(
    filename: str,
    stream,
    num_batches: Optional[int] = None,
    legacy: bool = False,
    resolution: TSRES = TSRES.NS1,
    display: bool = True,
    allow_duplicates: bool = True,
):
    """Block streaming variant of 'write_a1'.

    Input stream assumed to have TSRES.NS1 resolution.

    Args:
        filename: Destination filename.
        stream: Input stream, in TSRES.NS1 resolution.
        num_batches: Number of batches in stream.
        legacy: If output timestamp format should be in legacy format.
        resolution: Resolution of timestamps in input stream.
        display: If progress bar should be displayed during write.

    Note:
        Settling on a fixed resolution allows for more consistent interface
        for filtering functions. Code runtime with filtering will take a hit,
        however - this is a future todo.
    """
    if display:
        stream = tqdm.tqdm(stream, total=num_batches)
    with open(filename, "wb") as f:
        for t, p in stream:
            events = _consolidate_events(
                t, p, resolution, allow_duplicates=allow_duplicates
            )
            for line in events:
                if legacy:
                    line = int(line)
                    line = ((line & 0xFFFFFFFF) << 32) + (line >> 32)
                f.write(struct.pack("=Q", line))


def swrite_a0(
    filename: str,
    stream,
    num_batches: Optional[int] = None,
    legacy: Optional[bool] = None,
    resolution: TSRES = TSRES.NS1,
    display: bool = True,
    allow_duplicates: bool = True,
):
    """See documentation for 'swrite_a1'."""
    if display:
        stream = tqdm.tqdm(stream, total=num_batches)
    with open(filename, "w") as f:
        for t, p in stream:
            events = _consolidate_events(
                t, p, resolution, allow_duplicates=allow_duplicates
            )
            data = np.empty((2 * events.size,), dtype=np.uint32)
            data[0::2] = events & 0xFFFFFFFF
            data[1::2] = events >> 32
            for line in data:
                f.write(f"{line:08x}\n")


def swrite_a2(
    filename: str,
    stream,
    num_batches: Optional[int] = None,
    legacy: Optional[bool] = None,
    resolution: TSRES = TSRES.NS1,
    display: bool = True,
    allow_duplicates: bool = True,
):
    """See documentation for 'swrite_a1'."""
    if display:
        stream = tqdm.tqdm(stream, total=num_batches)
    with open(filename, "w") as f:
        for t, p in stream:
            data = _consolidate_events(
                t, p, resolution, allow_duplicates=allow_duplicates
            )
            for line in data:
                f.write(f"{line:016x}\n")


##############
#  PRINTERS  #
##############


def print_statistics(filename: str, t: list, p: list):
    """Prints statistics using timestamp event readers.

    Note:
        Maintained only for legacy reasons to support older
        reading mechanisms, e.g. 'read_a0'.
    """
    # Collect statistics
    count = np.count_nonzero
    _ch1 = (count(p & 0b0001 != 0),)
    _ch2 = (count(p & 0b0010 != 0),)
    _ch3 = (count(p & 0b0100 != 0),)
    _ch4 = (count(p & 0b1000 != 0),)
    print_statistics_report(
        filename=filename,
        num_events=len(t),
        num_detections=_ch1 + _ch2 + _ch3 + _ch4,
        ch1_counts=_ch1,
        ch2_counts=_ch2,
        ch3_counts=_ch3,
        ch4_counts=_ch4,
        multi_counts=count(np.isin(p, (0, 1, 2, 4, 8), invert=True)),
        non_counts=count(p == 0),
        start_timestamp=t[0],
        end_timestamp=t[-1],
        patterns=sorted(np.unique(p)),
    )


def print_statistics_stream(
    filename: str,
    stream,
    num_batches: Optional[int] = None,
    resolution: TSRES = TSRES.NS1,
    display: bool = True,
):
    """Prints statistics using timestamp event streamers."""
    # Calculate statistics
    first_t = None
    last_t = None
    num_events = 0
    num_detections = 0
    count_p1 = 0
    count_p2 = 0
    count_p3 = 0
    count_p4 = 0
    count_mp = 0
    count_np = 0
    set_p = set()

    # Processs batches of events
    count = np.count_nonzero
    if display:
        stream = tqdm.tqdm(stream, total=num_batches)
    for t, p in stream:
        _ch1 = count(p & 0b0001 != 0)
        _ch2 = count(p & 0b0010 != 0)
        _ch3 = count(p & 0b0100 != 0)
        _ch4 = count(p & 0b1000 != 0)
        num_events += len(p)
        num_detections += _ch1 + _ch2 + _ch3 + _ch4
        count_p1 += _ch1
        count_p2 += _ch2
        count_p3 += _ch3
        count_p4 += _ch4
        count_mp += count(np.isin(p, (0, 1, 2, 4, 8), invert=True))
        count_np += count(p == 0)
        set_p.update(np.unique(p))

        # Store timing only if timestamps not filtered out
        if len(t) != 0:
            if first_t is None:
                first_t = t[0] / resolution.value  # convert to nanoseconds
            last_t = t[-1] / resolution.value

    print_statistics_report(
        filename=filename,
        num_events=num_events,
        num_detections=num_detections,
        ch1_counts=count_p1,
        ch2_counts=count_p2,
        ch3_counts=count_p3,
        ch4_counts=count_p4,
        multi_counts=count_mp,
        non_counts=count_np,
        start_timestamp=first_t,
        end_timestamp=last_t,
        patterns=sorted(set_p),
    )


def print_statistics_report(
    filename: str,
    num_events: int,
    num_detections: int,
    ch1_counts: Optional[int] = None,
    ch2_counts: Optional[int] = None,
    ch3_counts: Optional[int] = None,
    ch4_counts: Optional[int] = None,
    multi_counts: Optional[int] = None,
    non_counts: Optional[int] = None,
    start_timestamp: Optional[float] = None,
    end_timestamp: Optional[float] = None,
    patterns: Optional[list] = None,
):
    """Prints the statistics report.

    All optional fields must be present if 'num_events' > 0.

    Events and detections have the following definitions: multi-channel events
    are counted as a single event for the former, and as separate events for
    the latter.
    """

    print(f"Name: {str(filename)}")
    if pathlib.Path(filename).is_file():
        filesize = pathlib.Path(filename).stat().st_size / (1 << 20)
        print(f"Filesize (MB): {filesize:.3f}")
    width = 0
    if num_events != 0:
        width = int(np.floor(np.log10(num_events))) + 1
    print(f"Total events    : {num_events:>{width}d}")
    if num_events != 0:
        duration = (end_timestamp - start_timestamp) * 1e-9
        print(f"  Detections    : {num_detections:>{width}d}")
        print(f"  Channel 1     : {ch1_counts:>{width}d}")
        print(f"  Channel 2     : {ch2_counts:>{width}d}")
        print(f"  Channel 3     : {ch3_counts:>{width}d}")
        print(f"  Channel 4     : {ch4_counts:>{width}d}")
        print(f"  Multi-channel : {multi_counts:>{width}d}")
        print(f"  No channel    : {non_counts:>{width}d}")
        print(f"Duration (s) : {duration:15.9f}")
        print(f"  ~ start    : {start_timestamp*1e-9:>15.9f}")
        print(f"  ~ end      : {end_timestamp*1e-9:>15.9f}")
        print(f"Event rate (/s)  : {int(num_events//duration)}")
        print(f"  Detection rate : {int(num_detections//duration)}")
        print(f"Detection patterns: {sorted(map(int, patterns))}")


###############
#  FILTERING  #
###############


def get_pattern_mask(
    p: list,
    pattern: int,
    mask: bool = False,
    invert: bool = False,
):
    """Returns a mask with bits set where patterns match, as well as result.

    The function behaves differently when the pattern is either used as
    a fixed pattern (when 'mask' is False), or as a bitmask (when 'mask'
    is True).

    When 'mask' is False, pattern matching with select events whose detector
    patterns match the exact pattern. Inverting with 'invert' as False will
    remove these events instead (i.e. pattern blacklist).

    When 'mask' is True, only events containing any bits in the pattern will
    be selected. Inverting with 'invert' as True will remove detector bits
    corresponding to any of said bits.

    Args:
        p: List of detector patterns.
        pattern: Pattern to match.
        mask: Whether pattern should be used as a bitmask.
        invert: Whether selection rules should be inverted. See documentation.

    Examples:

        >>> p = np.array([0,1,2,3,4,5,6,7,8]); t = np.array(p)

        >>> mask1, p1 = get_pattern_mask(p, 0b0101, False, False)
        >>> list(t[mask1]) == list(p1) and list(p1) == [5]
        True

        >>> mask2, p2 = get_pattern_mask(p, 0b0101, False, True)
        >>> list(t[mask2]) == list(p2) and list(p2) == [0,1,2,3,4,6,7,8]
        True

        >>> mask3, p3 = get_pattern_mask(p, 0b0101, True, False)
        >>> list(t[mask3]) == [1,3,4,5,6,7] and list(p3) == [1,1,4,5,4,5]
        True

        >>> mask4, p4 = get_pattern_mask(p, 0b0101, True, True)
        >>> list(t[mask4]) == [0,2,3,6,7,8] and list(p4) == [0,2,2,2,2,8]
        True

    Note:
        The function is designed to return a bitmask instead of directly
        filtering the timestamps, to allow subsequent post-processing based
        on the bitmask, e.g. inverting the returned bitmask, or mark patterns.
    """

    # Fixed pattern
    if not mask:
        pmask = p == pattern
        if invert:
            pmask = ~pmask
        masked = p[pmask]

    # Pattern as bitmask
    elif not invert:
        _p = p & pattern  # patterns containing bitmask bits
        pmask = _p != 0
        masked = _p[pmask]
    else:
        # Set bit 4 to indicate dummy events, since
        # non-events already do not contain the bit pattern
        _p = np.where(p == 0, 16, p)
        _p = _p ^ (_p & pattern)  # patterns with bitmask removed
        pmask = _p != 0
        masked = _p[pmask]  # return detector patterns after masking
        masked = np.where(masked == 16, 0, masked)  # recover dummy events

    return pmask, masked


def get_timing_bounds(
    t: list,
    start: Optional[float] = None,
    end: Optional[float] = None,
    resolution: TSRES = TSRES.NS1,
) -> list:
    """Returns indices of start and end for region where timings match.

    The timing array is already assumed to be sorted, as part of the timestamp
    filespec. If 'start' or 'end' is None, the range is assumed unbounded in respective
    direction. Start and end time are in seconds (not ns since typical usecase as rough
    filter anyway).

    This function is preferred over 'get_timing_mask' inline filtering,
    due to the overheads of masking.

    Args:
        t: Timestamp array, in units of 'resolution'.
        start: Start timestamp, in seconds.
        end: End timestamp, in seconds.
        resolution: Resolution of timestamps in timestamp array.

    Note:
        # Masking overhead
        >>> timeit
        >>> setup = "import numpy as np; size=1_000_000; a = np.random.random(size);"
        >>> setup_w_mask = setup + "mask = np.zeros(size).astype(bool);"
        >>> setup_w_mask += "mask[1000:999000] = True;"
        >>> timeit.timeit('a[1000:999000]', setup, number=1000)
        0.00011489982716739178
        >>> timeit.timeit('a[mask]', setup_w_mask, number=1000)
        1.7692412599921226
    """
    # Nothing in array
    NULL_RESULT = (0, 0)
    if len(t) == 0:
        return NULL_RESULT

    if start is not None:
        start *= 1e9 * resolution.value  # convert to same base for comp.
    if end is not None:
        end *= 1e9 * resolution.value
    _start = t[0]
    _end = t[-1]

    # Find start index
    if (start is None) or (start <= _start):
        i = 0
    elif start > _end:
        return NULL_RESULT
    else:
        i = bisect.bisect_left(t, start)  # binary search

    # Find end index
    if (end is None) or (end >= _end):
        j = len(t)
    elif end < _start:
        return NULL_RESULT
    else:
        j = bisect.bisect_right(t, end)

    return i, j


def get_timing_mask(
    t: list,
    start: Optional[float] = None,
    end: Optional[float] = None,
    resolution: TSRES = TSRES.NS1,
) -> list:
    """Returns a mask where timestamps are bounded between start and end.

    The timing array is already assumed to be sorted, as part of the timestamp
    filespec. If 'start' or 'end' is None, the range is assumed unbounded in respective
    direction. Start and end time are in seconds (not ns since typical usecase as rough
    filter anyway).

    Args:
        t: Timestamp array, in units of 'resolution'.
        start: Start timestamp, in seconds.
        end: End timestamp, in seconds.
        resolution: Resolution of timestamps in timestamp array.

    Examples:

        >>> to_ns = lambda ts: [t * 1e9 for t in ts]
        >>> t = np.array(to_ns([2,4,5,8]))

        >>> mask1 = get_timing_mask(t, 4, 8)  # range at boundaries
        >>> list(t[mask1]) == to_ns([4,5,8])
        True

        >>> mask2 = get_timing_mask(t, 3, 6)  # range outside boundaries
        >>> list(t[mask2]) == to_ns([4,5])
        True

        >>> mask3 = get_timing_mask(t, 9, 100)  # too high
        >>> list(t[mask3]) == to_ns([])
        True

        >>> mask4 = get_timing_mask(t, 0, 1)  # too low
        >>> list(t[mask4]) == to_ns([])
        True

        >>> mask5 = get_timing_mask(t, 5, 5)  # end is inclusive
        >>> list(t[mask5]) == to_ns([5])
        True

        >>> mask6 = get_timing_mask(t, start=3, end=None)  # only lower bound
        >>> list(t[mask6]) == to_ns([4,5,8])
        True

        >>> mask7 = get_timing_mask(t, start=None, end=4)  # only upper bound
        >>> list(t[mask7]) == to_ns([2,4])
        True
    """
    mask = np.zeros(len(t), dtype=bool)
    i, j = get_timing_bounds(t, start, end, resolution)
    mask[i:j] = True
    return mask


##########
#  MAIN  #
##########


def generate_parser():
    parser = argparse.ArgumentParser(
        description="Converts between different timestamp7 formats"
    )

    parser.add_argument(
        "-q", "--quiet", action="store_true", help="Suppress progress indicators"
    )
    # Support for older read-write mechanisms, i.e. disable batch streaming
    parser.add_argument(
        "--inmemory",
        action="store_true",
        help="Disable batch streaming (retained for legacy reasons)",
    )
    parser.add_argument(
        "--rollover",
        action="store_true",
        help="Include rollover dummy events (generally unwanted)",
    )

    pgroup = parser.add_argument_group("Timestamp formats")
    pgroup.add_argument(
        "-A",
        choices=["0", "1", "2"],
        default="1",
        help="Input timestamp format (default: 1)",
    )
    pgroup.add_argument(
        "-X", action="store_true", help="Use legacy input format (default: False)"
    )
    pgroup.add_argument(
        "-a",
        choices=["0", "1", "2"],
        default="1",
        help="Output timestamp format (default: 1)",
    )
    pgroup.add_argument(
        "-x", action="store_true", help="Use legacy output format (default: False)"
    )

    # Filtering
    pgroup = parser.add_argument_group("Pattern filtering")
    pgroup.add_argument(
        "--pattern",
        type=int,
        metavar="",
        help="Specify pattern for filtering (e.g. 5 for ch1+ch3)",
    )
    pgroup.add_argument(
        "--mask",
        action="store_true",
        help="Use pattern as mask instead of fixed (default: False)",
    )
    pgroup.add_argument(
        "--invert",
        action="store_true",
        help="Exclude pattern instead of include (default: False)",
    )

    pgroup = parser.add_argument_group("Timestamp filtering")
    pgroup.add_argument(
        "--start", type=float, metavar="", help="Specify start timestamp, in seconds"
    )
    pgroup.add_argument(
        "--end", type=float, metavar="", help="Specify end timestamp, in seconds"
    )
    pgroup.add_argument(
        "--relative-time",
        action="store_true",
        help="Use relative time to start timestamp instead.",
    )

    parser.add_argument("infile", help="Input timestamp file")
    parser.add_argument(
        "outfile",
        nargs="?",
        const="",
        help="Output timestamp file (optional: not required for printing)",
    )
    return parser


def main():  # noqa: PLR0915
    # Print help if no arguments supplied
    parser = generate_parser()
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    args = parser.parse_args()
    print_only = args.outfile is None
    ignore_rollover = not args.rollover

    read = [read_a0, read_a1, read_a2][int(args.A)]
    sread = [sread_a0, sread_a1, sread_a2][int(args.A)]
    write = [write_a0, write_a1, write_a2][int(args.a)]
    swrite = [swrite_a0, swrite_a1, swrite_a2][int(args.a)]

    # Check file size
    filepath = pathlib.Path(args.infile)
    if not filepath.is_file():
        raise ValueError(f"'{args.infile}' is not a file.")

    # Define filter function for event and event stream
    has_filtering = (
        (args.pattern is not None) or (args.start is not None) or (args.end is not None)
    )

    # Apply relative time
    if args.relative_time and args.end is not None:
        first_t, _ = read_a1_start_end(
            filepath,
            args.X,
            ignore_rollover=ignore_rollover,
        )
        args.end += first_t  # default 1ns resolution

    def event_filter(t, p, resolution: TSRES = TSRES.NS1):
        if args.pattern is not None:
            mask, p = get_pattern_mask(p, args.pattern, args.mask, args.invert)
            t = t[mask]
        if args.start is not None or args.end is not None:
            mask = get_timing_mask(t, args.start, args.end, resolution)
            t = t[mask]
            p = p[mask]
        return t, p

    def inline_filter(stream, resolution: TSRES = TSRES.NS1):
        commenced_reading = False
        for t, p in stream:
            _t, _p = event_filter(t, p, resolution)
            if commenced_reading and len(_t) == 0:
                break
            elif len(_t) != 0:
                commenced_reading = True
            yield _t, _p

    # Use legacy read-write mechanisms
    if args.inmemory:
        t, p = read(filepath, args.X, ignore_rollover=ignore_rollover)
        t, p = event_filter(t, p)
        if print_only:
            print_statistics(filepath, t, p)
        else:
            write(args.outfile, t, p, args.x)

    # Check if printing stream first
    elif print_only:
        if has_filtering:
            stream, num_batches = sread(
                filepath, args.X, ignore_rollover=ignore_rollover
            )  # default 1ns floating-point resolution
            stream = inline_filter(stream)
            print_statistics_stream(
                filepath, stream, num_batches, display=(not args.quiet)
            )
        else:
            # Disable timestamp formatting to speed up reads
            stream, num_batches = sread(
                filepath, args.X, TSRES.PS4, False, ignore_rollover=ignore_rollover
            )
            print_statistics_stream(
                filepath,
                stream,
                num_batches,
                resolution=TSRES.PS4,
                display=(not args.quiet),
            )

    # Write out
    else:
        stream, num_batches = sread(
            filepath, args.X, TSRES.PS4, False, ignore_rollover=ignore_rollover
        )
        stream = inline_filter(stream, TSRES.PS4)
        swrite(
            args.outfile,
            stream,
            num_batches,
            args.x,
            resolution=TSRES.PS4,
            display=(not args.quiet),
        )


if __name__ == "__main__":
    main()
