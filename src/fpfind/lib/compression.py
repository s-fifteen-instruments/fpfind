#!/usr/bin/env python3
"""Stores routines relevant for timestamp data compression and archival.

Performing delta-encoding can shrink compatible timestamp datasets by
approximately 50%, and takes a relatively fast amount of time compared
to compressing directly using DEFLATE/lzma, i.e. around 6s per GB of
timestamp data on AMD 5500U. Internally uses numpy's '.npz' file
specification to store the timestamps directly as integer arrays into
the file.

This delta-encoded timestamps can be compressed further if one needs to
save even more space, up to x1/3 of original size. See examples below for
performance comparisions.

Examples:

    >>> import lzma
    >>> import pathlib
    >>> from fpfind.lib.compression import DTS
    >>> def compress_xz(filename):
    ...     target = filename + ".xz"
    ...     with open(filename, "rb") as fin, lzma.open(fout, "w") as lzf:
    ...         while (buffer := fin.read(1_000_000)) != b"":
    ...             lzf.write(buffer)
    >>> def get_sizes(filename):
    ...     size = lambda path: pathlib.Path(path).stat().st_size
    ...     print(" raw:", size(filename))
    ...     print("gzip:", size(filename + ".gz"))
    ...     print("  xz:",size(filename + ".xz"))

    # Original timestamp data
    >>> compress_xz("timestamps.a1.ts")
    >>> get_sizes("timestamps.a1.ts")
     raw: 100000000
    gzip: 54689009
      xz: 42081640

    # Delta-encoded timestamp data
    >>> DTS.encode("timestamps.a1.ts")
    >>> compress_xz("timestamps.a1.dts")
    >>> get_sizes("timestamps.a1.dts")
     raw: 50000526
    gzip: 40258198
      xz: 33795448

    # Load from encoded data
    >>> ts, ps = DTS.load("timestamps.a1.dts")

TODO:
    Add support for streaming archival for large datasets (that cannot
    fit into memory).
"""

import dataclasses
import pathlib
from typing import TYPE_CHECKING, Optional, Tuple

import numpy as np
import numpy.typing as npt
from numpy.typing import NDArray

import fpfind.lib.parse_timestamps as parser
from fpfind import TSRES
from fpfind.lib.typing import (
    DetectorArray,
    Integer,
    PathLike,
    TimestampArray,
    TimestampDetectorStream,
)

DECODED_FILE_EXT: str = ".ts"
ENCODED_FILE_EXT: str = ".dts"


@dataclasses.dataclass
class DTS:
    """Wrapper for delta-encoded timestamps.

    Timestamps can be stored raw, where 'data' represents 32-bit deltas and 'ps'
    represents the 6-bit timestamp metadata (including blinding/dummy), or stored
    compressed, where 'data' stores 28-bit deltas and 4-bit detector pattern and
    'ps' is uninitialized. Use the 'lossless' parameter to disable compression.

    Examples:
        >>> DTS.encode("mytimestamps.a1.ts")
        >>> ts, ps = DTS.load("mytimestamps.a1.dts", fractional=False)
    """

    ts0: np.uint64
    const: np.uint64
    data: npt.NDArray[np.uint32]
    ps: Optional[npt.NDArray[np.uint8]] = None

    @staticmethod
    def encode(
        filename: PathLike,
        legacy: bool = False,
        lossless: bool = False,
        target: Optional[PathLike] = None,
    ) -> None:
        """Encodes a timestamp dataset into a delta-encoding variant.

        If 'target' is not specified, the encoded file will be co-located with the
        timestamp dataset with the '.dts' file extension.

        Args:
            filename: Path to timestamp dataset.
            legacy: Whether timestamp has legacy encoding.
            lossless: Whether to use lossless encoding at the expense of space.
            target: Write path to encoded timestamp dataset.
        """
        if target is None:
            target = pathlib.Path(filename).with_suffix(ENCODED_FILE_EXT)
        data = DTS._load_decoded(filename, legacy=legacy, lossless=lossless)
        if not lossless:
            data._compress()
        data._save_encoded(target)

    @staticmethod
    def decode(
        filename: PathLike, legacy: bool = False, target: Optional[PathLike] = None
    ) -> None:
        """Decodes a delta-encoded timestamp dataset into the original dataset.

        Basically the inverse of 'DTS.encode()'.

        This is not a 1-to-1 decoding (even more so if lossless encoding was not used).
        This can be implemented by reinserting the constant field, as well as preserving
        the dummy/blinding flags.
        """
        if target is None:
            target = pathlib.Path(filename).with_suffix(DECODED_FILE_EXT)
        data = DTS._load_encoded(filename)
        data._save_decoded(target, legacy=legacy)

    @staticmethod
    def load(
        filename: PathLike, resolution: TSRES = TSRES.NS1, fractional: bool = True
    ) -> Tuple[TimestampArray, DetectorArray]:
        """Convenience function to load timestamps from encoded dataset.

        Functions similar to 'fpfind.lib.parse_timestamps.read_a1', read the
        corresponding documentation for more details on the parameters.

        Args:
            filename: Path to encoded timestamps.
            resolution: Target resolution of decoded timestamps.
            fractional: Whether to preserve fractional parts.
        """
        data = DTS._load_encoded(filename)
        ts, ps = data._get(raw=False)
        if TYPE_CHECKING:  # for type casting only
            ts = np.asarray(ts, dtype=np.uint64)
        ts = parser._format_timestamps(ts, resolution=resolution, fractional=fractional)
        return ts, ps

    ##############
    #  INTERNAL  #
    ##############

    def _get(self, raw: bool = True) -> Tuple[TimestampArray, DetectorArray]:
        """Returns raw ts and ps values."""
        dts = self.data
        ps = self.ps
        if ps is None:
            ps = self.data & np.uint32(0xF)
            dts = self.data >> np.uint32(4)
        elif not raw:
            ps = ps & np.uint8(0xF)

        ts = np.cumsum(dts, dtype=np.uint64) + self.ts0
        ps = np.asarray(ps, dtype=np.uint32)
        return ts, ps

    def _compress(self) -> None:
        """Perform in-place lossy compression.

        Using 32-bit constructs similar to readevents7 fast mode '-f', but with larger
        timestamp width. Consists of 28-bit timing and 4-bit pattern, which corresponds
        to 2^20ns (~1ms) for 4ps timing resolution.

        This is hidden from user to avoid dummy events leaking into dataset, when
        compression is performed on dataset read with 'ignore_rollover=False'.

        Note:
            Future implementations can consider dynamically storing events for
            preserving rare instances of large timing differences.
        """
        if self.ps is None:
            return  # already in lossy format

        # Compress dts into 28-bit widths
        dts = self.data
        dts_mask = (1 << 28) - 1
        if not np.all(dts <= dts_mask):
            raise ValueError("Timestamp deltas longer than 1ms are not supported.")
        dts &= np.uint32(dts_mask)

        # Remove blinding / dummy event information
        ps = self.ps
        ps &= np.uint8(0xF)
        events = ((dts << np.uint32(4)) + ps).astype(np.uint32)

        # Modify in-place
        self.data = events
        self.ps = None

    def _save_encoded(self, filename: PathLike) -> None:
        """
        Can consider variable width in blocks of 8-bits, i.e. storing data with dtype
        of np.uint8, then expand accordingly to save marginal amounts of space
        (corresponding to 1/64 of original filesize per bit).
        """
        metadata = np.array([self.ts0, self.const], dtype=np.uint64)
        kwargs = {}
        if self.ps is not None:
            kwargs = {"ps": self.ps}

        target = pathlib.Path(filename).with_suffix(".npz")
        np.savez(
            target,  # '.npz' suffix is forced by numpy
            allow_pickle=False,
            data=self.data,
            metadata=metadata,
            **kwargs,
        )
        target.rename(filename)  # rename

    @staticmethod
    def _load_encoded(filename: PathLike) -> "DTS":
        with np.load(filename) as d:
            ts0, const = d["metadata"]
            return DTS(ts0, const, d["data"], d.get("ps", None))

    def _save_decoded(self, filename: PathLike, legacy: bool = False) -> None:
        ts, ps = self._get()
        parser.write_a1(filename, ts, ps, legacy=legacy, resolution=parser.TSRES.PS4)

    @staticmethod
    def _load_decoded(
        filename: PathLike,
        legacy: bool = False,
        ignore_const: bool = True,
        lossless: bool = False,
    ) -> "DTS":
        """Constant field is automatically ignored."""
        ignore_rollover = not lossless
        ts, ps = _read_a1_lossless(filename, legacy, ignore_rollover=ignore_rollover)
        ts = np.asarray(ts, dtype=np.uint64)
        ts0 = ts[0]

        # Verify constant field
        const = 0b1100
        if not ignore_const:
            consts = np.unique((ps & np.uint32(0x3C0)) >> np.uint32(6))
            assert len(consts) == 1
            const = consts[0]
            assert const == 0b1100
        const = np.uint64(const)

        ps = np.asarray(ps & np.uint32(0x3F), dtype=np.uint8)

        # Perform run-length encoding
        dts = np.diff(ts, prepend=[ts0])
        dts_mask = (1 << 32) - 1
        if not np.all(dts <= dts_mask):
            raise ValueError("Timestamp deltas longer than 17ms are not supported.")
        dts = dts.astype(np.uint32)  # 2^32 x 4ps = ~17ms

        data = DTS(ts0, const, dts, ps)
        if not lossless:
            data._compress()
        return data


def _read_a1_lossless(
    filename: PathLike,
    legacy: bool = False,
    resolution: TSRES = TSRES.PS4,
    fractional: bool = False,
    buffer_size: Integer = 100_000,
    ignore_rollover: bool = True,
) -> Tuple[TimestampArray, DetectorArray]:
    """
    Lossless variant of 'read_a1', where detector pattern is the full 10-bit
    timestamp metadata. Used mainly for binary compression, so most of the
    options carried over from 'read_a1' are fixed, including:

      * resolution = TSRES.PS4  # default to widest timestamp field
      * fractional = False      # retrieve exact integer deltas
      * buffer_size = 100_000   # in-memory reads are disabled for efficiency
      * ignore_rollover = True  # ignore unnecessary timestamps

    Precluded from main implementation to avoid bloat.
    """
    streamer, _ = _sread_a1_lossless(
        filename, legacy, resolution, fractional, buffer_size, ignore_rollover
    )
    ts = []
    ps = []
    commenced_reading = False
    for t, p in streamer:
        # Terminate if no more timings readable
        if commenced_reading and len(t) == 0:
            break
        elif len(t) != 0:
            commenced_reading = True
        ts.append(t)
        ps.append(p)

    # Return array of correct dtype if no data
    if len(ts) == 0:
        empty = np.array([], dtype=np.uint64)
        t = parser._format_timestamps(empty, resolution, fractional)
        p = np.array([], dtype=np.uint32)  # consistent with %I format code
        return t, p

    t = np.concatenate(ts)
    p = np.concatenate(ps)
    return t, p


def _sread_a1_lossless(
    filename: PathLike,
    legacy: bool = False,
    resolution: TSRES = TSRES.NS1,
    fractional: bool = True,
    buffer_size: Integer = 100_000,
    ignore_rollover: bool = False,
) -> Tuple[TimestampDetectorStream, int]:
    def _sread_a1():
        with open(filename, "rb") as f:
            while True:
                buffer = f.read(buffer_size * 8)  # 8 bytes per event
                if len(buffer) == 0:
                    break
                data = np.frombuffer(buffer, dtype="<u4").reshape(-1, 2)
                yield _parse_a1_lossless(
                    data, legacy, resolution, fractional, ignore_rollover
                )

    size = pathlib.Path(filename).stat().st_size
    num_batches = int(((size - 1) // (buffer_size * 8)) + 1)
    return _sread_a1(), num_batches


def _parse_a1_lossless(
    data: NDArray[np.uint32],
    legacy: bool = False,
    resolution: TSRES = TSRES.NS1,
    fractional: bool = True,
    ignore_rollover: bool = False,
) -> Tuple[TimestampArray, DetectorArray]:
    assert data.ndim == 2 and data.shape[1] == 2

    high_pos = 1
    low_pos = 0
    if legacy:
        high_pos, low_pos = low_pos, high_pos

    if ignore_rollover:
        r = (data[:, low_pos] & 0b10000).astype(bool)  # check rollover flag
        data = data[~r]

    high_words = data[:, high_pos].astype(np.uint64) << np.uint(22)
    low_words = data[:, low_pos] >> np.uint(10)
    ts = parser._format_timestamps(high_words + low_words, resolution, fractional)
    ps = data[:, low_pos] & np.uint32(0x3FF)  # only this field is changed
    return ts, ps
