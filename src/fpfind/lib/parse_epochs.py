#!/usr/bin/env python3
"""Takes in an epoch and parses it into timestamps.

Currently only supports T1 and T2 epoch unpacking.

References:
    [1] File specification: https://github.com/s-fifteen-instruments/qcrypto/blob/Documentation/docs/source/file%20specification.rst
    [2] Epoch header definitions: https://github.com/s-fifteen-instruments/QKDServer/blob/master/S15qkd/utils.py

Changelog:
    2023-02-15 Justin: Init
"""

import datetime as dt
import pathlib
from dataclasses import dataclass, field
from struct import pack, unpack
from typing import NamedTuple

import numpy as np

from fpfind import NP_PRECISEFLOAT, TSRES
from fpfind.lib._logging import get_logger

logger = get_logger(__name__, level="warning")


class HeadT1(NamedTuple):
    tag: int
    epoch: int
    length_bits: int
    bits_per_entry: int
    base_bits: int


class HeadT2(NamedTuple):
    tag: int
    epoch: int
    length_bits: int
    timeorder: int
    base_bits: int
    protocol: int


class HeadT3(NamedTuple):
    tag: int
    epoch: int
    length_entry: int
    bits_per_entry: int


class HeadT4(NamedTuple):
    tag: int
    epoch: int
    length: int
    timeorder: int
    base_bits: int


# Header readers
# Extracted from ref [2]


def read_T1_header(file_name: str):
    file_name = str(file_name)
    with open(file_name, "rb") as f:
        head_info = f.read(4 * 5)
    headt1 = HeadT1._make(unpack("iIIii", head_info))
    if headt1.tag not in (0x101, 1):
        logger.warning(f"{file_name} is not a Type2 header file")
    if hex(headt1.epoch) != ("0x" + file_name.split("/")[-1]):
        logger.warning(
            f"Epoch in header {headt1.epoch} does not match epoc filename {file_name}"
        )
    return headt1


def read_T2_header(file_name: str):
    file_name = str(file_name)
    with open(file_name, "rb") as f:
        head_info = f.read(4 * 6)
    headt2 = HeadT2._make(unpack("iIIiii", head_info))
    if headt2.tag not in (0x102, 2):
        logger.warning(f"{file_name} is not a Type2 header file")
    if hex(headt2.epoch) != ("0x" + file_name.split("/")[-1]):
        logger.warning(
            f"Epoch in header {headt2.epoch} does not match epoc filename {file_name}"
        )
    return headt2


def read_T3_header(file_name: str):
    file_name = str(file_name)
    with open(file_name, "rb") as f:
        head_info = f.read(4 * 4)
    headt3 = HeadT3._make(unpack("iIIi", head_info))
    if headt3.tag not in (0x103, 3):
        logger.warning(f"{file_name} is not a Type3 header file")
    if hex(headt3.epoch) != ("0x" + file_name.split("/")[-1]):
        logger.warning(
            f"Epoch in header {headt3.epoch} does not match epoc filename {file_name}"
        )
    return headt3


def read_T4_header(file_name: str):
    file_name = str(file_name)
    with open(file_name, "rb") as f:
        head_info = f.read(4 * 5)
    headt4 = HeadT4._make(unpack("IIIII", head_info))
    if headt4.tag not in (0x104, 4):
        logger.warning(f"{file_name} is not a Type4 header file")
    if hex(headt4.epoch) != ("0x" + file_name.split("/")[-1]):
        logger.warning(
            f"Epoch in header {headt4.epoch} does not match epoc filename {file_name}"
        )
    return headt4


# Epoch utility functions


def date2epoch(datetime=None):
    """Returns current epoch number as hex.

    Epochs are in units of 2^32 * 125e-12 == 0.53687 seconds.

    Args:
        datetime: Optional datetime to convert to epoch.
    """
    if datetime is None:
        datetime = dt.datetime.now()

    total_seconds = int(datetime.timestamp())
    epoch_val = int(total_seconds / 125e-12) >> 32
    return f"{epoch_val:08x}"


def epoch2date(epoch):
    """Returns datetime object corresponding to epoch.

    Epoch must be in hex format, e.g. "ba1f36c0".

    Args:
        epoch: Epoch in hex format, e.g. "ba1f36c0".
    """
    total_seconds = (int(epoch, base=16) << 32) * 125e-12
    return dt.datetime.fromtimestamp(total_seconds)


def epoch2int(epoch):
    return int(epoch, base=16)


def int2epoch(value):
    return hex(value)[2:]


# Epoch readers
# Implemented as per filespec [1]


def read_T1(
    filename: str,
    full_epoch: bool = False,
    resolution: TSRES = TSRES.NS1,
    fractional: bool = True,
):
    """Returns timestamp and detector information from T1 files.

    Timestamp and detector information follow the formats specified in
    'parse_timestamps' module, for interoperability.

    By default, the function returns timestamps in units of 1ns with as
    high a resolution as possible (i.e. including fractional values).

    TODO

    'full' is False by default, which indicates only the raw-timestamp
    -equivalent information is returned, i.e. timestamp in units of 1ns
    stored as a 64-bit integer.
    'full' set to True returns full timestamps (timestamps + full epoch
    information) stored as 128-bit float.
    The reason for this difference is due to the full timestamp taking up
    more than 64-bits equivalent information in units of 1ns, so only a
    float can support this dynamic range.

    Note 128-bit float has precision of 113 bits, while 64-bit float only has
    53-bit precision (insufficient resolution).

    If 'raw' is True, timestamps are stored in units of 4ps instead of 1ns.

    Raises:
        ValueError: Number of epochs does not match length in header.
    """
    # Header details
    header = read_T1_header(filename)
    length = header.length_bits

    timestamps = []
    bases = []

    # For each timestamp event, 54-bits correspond to timestamp,
    # of which 17-bit MSB is LSB of epoch, 37-bit LSB is 125ps resolution timestamp
    with open(filename, "rb") as f:
        f.read(5 * 4)  # dump header
        while True:
            event_high = f.read(4)
            event_low = f.read(4)

            value = (
                int.from_bytes(event_high, byteorder="little") << 32
            ) + int.from_bytes(event_low, byteorder="little")

            if value == 0:
                break  # termination

            # Timestamp in units of 4ps, stored in 54-bit MSB
            timestamp = value >> 10
            timestamps.append(timestamp)

            base = value & 0b1111
            bases.append(base)

    # Validity checks
    if length and len(timestamps) != length:
        raise ValueError("Number of epochs do not match length specified in header")

    # Add epoch if required, at most 32-17 = 15 bits
    epoch_msb = header.epoch >> 17

    # Right now in units of 4ps
    # Check if need to store floating point
    if fractional:
        # Convert to desired resolution
        timestamps = np.array(timestamps, dtype=NP_PRECISEFLOAT)
        timestamps = timestamps / (TSRES.PS4.value / resolution.value)
        epoch_msb = NP_PRECISEFLOAT(epoch_msb << 54)
        epoch_msb = epoch_msb / (TSRES.PS4.value / resolution.value)

    # Python integer objects can be arbitrarily long
    elif resolution in (TSRES.PS4,) and full_epoch:
        timestamps = np.array(timestamps, dtype=object)
        epoch_msb = int(epoch_msb) << 54

    # Everything of resolution 125 ps and larger can fit
    # in 64-bit unsigned integer, i.e. PS125, NS1.
    else:
        bitdiff = round(np.log2(TSRES.PS4.value / resolution.value))
        timestamps = np.array(timestamps, dtype=np.uint64)
        timestamps = timestamps >> bitdiff
        epoch_msb = np.uint64(epoch_msb) << np.uint64(54 - bitdiff)

    # Add epoch
    if full_epoch:
        timestamps += epoch_msb

    return timestamps, np.array(bases)


def extract_bits(msb_size: int, buffer: int, size: int, fileobject=None):
    """Return 'msb_size'-bits from MSB of buffer, and the rest.

    If insufficient bits to load value, data is automatically loaded
    from file 'fileobject' provided (must be already open in "rb" mode).

    Raises:
        ValueError: msb_size > size but fileobject not provided.

    Example:

        >>> msb_size = 4
        >>> buffer = 0b0101110; size = 7
        >>> extract_bits(msb_size, buffer, size)
        (5, 6, 3)

        Explanation:
          - msb => 0b0101 = 4
          - buffer => 0b110 = 5

    Note:
        Integrating with 'append_word_from_file' since this
        operation always precedes the value retrieval. This function
        can be reused in bit-packing schemes of other epoch types.
    """
    # Append word from fileobject if insufficient bits
    if size < msb_size and fileobject:
        buffer <<= 32
        buffer += int.from_bytes(fileobject.read(4), byteorder="little")
        # buffer += fileobject.get_next_int()
        size += 32

    # Extract MSB from buffer
    msb = buffer >> (size - msb_size)
    buffer &= (1 << (size - msb_size)) - 1
    size -= msb_size
    # logger.debug("Extracted %d bytes: %d", msb_size, msb)
    return msb, buffer, size


class BufferedFileRead:
    """Stores internal buffer to minimize overall disk reads.

    Deprecated.

    Note:
        Written to stand in for the 'fileobject' object in 'extract_bits', so
        that additional block-based buffering can be performed in-situ.
    """

    def __init__(self, f, buffer_size: int = 100_000):
        self.f = f
        self.buffer_size = buffer_size

        # Track intermediate buffer
        self.buffer = b""
        self.pos = 0

    def read(self, num_bytes: int):
        # Read desired amount of bytes
        result = self.buffer[self.pos : self.pos + num_bytes]
        self.pos += num_bytes
        if len(result) == num_bytes:
            return result

        # Insufficient data, load more data into buffer read from file if available
        self.buffer = self.f.read(self.buffer_size * 8)
        if len(self.buffer) == 0:
            return None

        # Note shortfall must be less than number of bytes loaded into buffer
        # TODO: Allow for more data than in buffer read size, i.e. while loop check.
        shortfall = num_bytes - len(result)
        result += self.buffer[:shortfall]
        self.pos = shortfall
        return result


class BufferedFileRead2:
    """Stores internal buffer to minimize overall disk reads.

    Deprecated.

    Here we use 'np.frombuffer' to bulk process bytestring reads into integers,
    then use a pointer to track the next integer to read. This is more similar
    to that in 'costream::get_stream_2'.

    Note:
        Written to stand in for the 'fileobject' object in 'extract_bits', so
        that additional block-based buffering can be performed in-situ.

    Usage:
        >>> f = BufferedFileRead(f)
        >>> f.get_next_int()  # replaces 'int.from_bytes(f.read(4), "little")'
    """

    def __init__(self, f, buffer_size: int = 100_000):
        self.f = f
        self.buffer_size = buffer_size

        # Track intermediate buffer
        self.buffer = []
        self.pos = 0

    def read(self, num_bytes) -> bytes:
        """Reads fixed number of bytes.

        Note:
            Do not use this function except for compatibility reasons - this
            performs a double conversion from integer then back to bytestring.
            This is done to align with the new internal buffering method, where
            bulk processing of 32-bit integers is performed by 'np.frombuffer'.

            This new method the costly call of 'int.from_bytes' on individual
            32-bits constructs. For an epoch file of 800k events, this improves
            the readtime from 3.7s to.

            Update: Fake news! Buffering doesn't really help much due to the
            additional overhead.
        """
        assert num_bytes % 4 == 0, "must be read in units of 32-bits"
        result = []
        for i in range(num_bytes // 4):
            result.append(self.get_next_int())
        result = b"".join(map(lambda v: v.to_bytes(4, "little"), [304, 23, 12, 2]))
        return result

    def get_next_int(self):
        # Insufficient data, load more data into buffer read from file if available
        if self.pos == len(self.buffer):
            self.buffer = np.frombuffer(self.f.read(self.buffer_size * 8), dtype="=I")
            self.pos = 0
            if len(self.buffer) == 0:
                return None

        result = self.buffer[self.pos]
        self.pos += 1
        return result


def read_T2(
    filename: str,
    full_epoch: bool = False,
    resolution: TSRES = TSRES.NS1,
    fractional: bool = True,
):
    """Returns timestamp and detector information from T2 files.

    Timestamp and detector information follow the formats specified in
    'parse_timestamps' module, for interoperability.

    Raises:
        ValueError: Length and termination bits inconsistent with filespec.
    """
    if resolution not in TSRES:
        raise ValueError("Timestamp resolution must be one of enumeration TSRES")

    # Header details
    header = read_T2_header(filename)
    timeorder = header.timeorder
    timeorder_extended = 32  # just being pedantic
    basebits = header.base_bits
    length = header.length_bits

    # Accumulators
    timestamps = []  # all timestamps, units of 125ps
    bases = []

    # Add 17-bit epoch LSB as per regular timestamp spec,
    # i.e. (17-bit epoch LSB)(32-bit timestamp) => 49-bits
    epoch_lsb = header.epoch & ((1 << 17) - 1)
    timestamp = np.uint64(epoch_lsb << 32)  # current timestamp, deltas stored in T2
    buffer = 0
    size = 0  # number of bits in buffer
    with open(filename, "rb") as f:
        # f = BufferedFileRead(f)
        f.read(24)  # remove header
        while True:
            # Read timing information
            timing, buffer, size = extract_bits(timeorder, buffer, size, f)

            # End-of-file check
            if timing == 1:
                base, buffer, size = extract_bits(basebits, buffer, size, f)
                if base != 0:
                    raise ValueError(
                        "File inconsistent with T2 epoch definition: "
                        "Terminal base bits not zero"
                    )
                break

            # Check if timing is actually in extended format, i.e. 32-bits
            if timing == 0:
                timing, buffer, size = extract_bits(timeorder_extended, buffer, size, f)

            timestamp += timing
            timestamps.append(timestamp)

            # Extract detector pattern, but not important here
            # Note we are guaranteed (timeorder+basebits < 32)
            base, buffer, size = extract_bits(basebits, buffer, size, f)
            bases.append(base)

    # Validity checks
    if length and len(timestamps) != length:
        raise ValueError("Number of epochs do not match length specified in header")

    # Add full epoch if required
    epoch_msb = header.epoch >> 17

    # Right now in units of 125ps
    # Check if need to store floating point
    if fractional:
        # Convert to desired resolution
        timestamps = np.array(timestamps, dtype=NP_PRECISEFLOAT)
        timestamps = timestamps / (TSRES.PS125.value / resolution.value)
        epoch_msb = NP_PRECISEFLOAT(epoch_msb << 49)
        epoch_msb = epoch_msb / (TSRES.PS125.value / resolution.value)

    # Python integer objects can be arbitrarily long
    # for resolutions smaller than 125ps, i.e. larger value
    elif resolution.value > TSRES.PS125.value:
        bitdiff = round(np.log2(resolution.value / TSRES.PS125.value))
        timestamps = np.array(timestamps, dtype=object)
        timestamps = timestamps << bitdiff
        epoch_msb = int(epoch_msb) << (49 + bitdiff)

    # Everything of resolution 125 ps and larger can fit
    # in 64-bit unsigned integer, i.e. PS125, NS1.
    else:
        bitdiff = round(np.log2(TSRES.PS125.value / resolution.value))
        timestamps = np.array(timestamps, dtype=np.uint64)
        timestamps = timestamps >> bitdiff
        epoch_msb = np.uint64(epoch_msb) << np.uint64(49 - bitdiff)

    # Add epoch
    if full_epoch:
        timestamps += epoch_msb

    return timestamps, np.array(bases)


def write_T1(directory, full_epoch, timestamps, detectors):
    # Fit epoch to within 32-bits LSB
    full_epoch &= (1 << 32) - 1

    directory = pathlib.Path(directory)
    target = directory / f"{full_epoch:x}"

    # Broadcast detectors if single value given
    if np.array(detectors).ndim == 0:
        detectors = np.ones(len(timestamps)) * detectors
    elif len(detectors) != len(timestamps):
        raise ValueError("Lengths of timestamps and detectors must match")

    with open(target, "wb") as f:
        header = pack("iIIii", 0x101, full_epoch, len(timestamps), 49, 4)
        f.write(header)

        # Write events, 4ps resolution
        # 17-bit LSB epoch fits into a timestamp event
        timestamp_epoch = (full_epoch & ((1 << 17) - 1)) << (32 + 5)
        for timestamp, detector in zip(timestamps, detectors):
            full_timestamp = timestamp_epoch + timestamp
            event = (full_timestamp << 10) + int(detector)

            # Write high word, then low word
            f.write(pack("<II", event >> 32, event & 0xFFFF_FFFF))

        # Termination
        f.write(pack("II", 0, 0))

    return target


@dataclass
class ServiceT3:
    """Dataclass for T3 file formats.

    Coincidence matrix mapped directly as 2D array:

    [
        (VV),  VA , (VH),  VD ,
         AV , (AA),  AH , (AD),
        (HV),  HA , (HH),  HD ,
         DV , (DA),  DH , (DD),
    ]
    """

    head: HeadT3
    coinc_matrix: list = field(default_factory=lambda: [0] * 16)
    garbage: list = field(default_factory=lambda: [0, 0])
    okcount: int = 0
    qber: float = 0.5  # default random


def read_T3(file_name: str) -> ServiceT3:
    """Used to process rawkey files.

    Example data:
    01000001_00100101_00100010_00010100
      \\    \\
     alice  bob
     0100   0001

    Note:
        Unpacking was done wrongly in original diagnosis.c code.

        Functionally the same as 'service_T3' in QKDServer.utils, with
        some syntatic comments.
    """
    # Map single-bit detections to index in 4-array for coinc matrix,
    # with multi-bit detections set to -1
    decode = [-1, 0, 1, -1, 2, -1, -1, -1, 3, -1, -1, -1, -1, -1, -1, -1]
    er_coinc_id = [0, 5, 10, 15]  # VV, AA, HH, DD
    gd_coinc_id = [2, 7, 8, 13]  # VH, AD, HV, DA
    detections = []  # store detection events of pairs
    headt3 = read_T3_header(file_name)
    service = ServiceT3(headt3)

    # Check if number of bits per entry according to filespec
    if headt3.bits_per_entry != 8:
        logger.warning("Not a service file with 8 bits per entry")

    # Digest file
    HEADER_INFO_SIZE = 16
    data_b = b"\x00\x00\x00\x00"
    with open(file_name, "rb") as f:
        f.seek(HEADER_INFO_SIZE)
        while True:
            word = f.read(4)
            if word == b"":
                break

            # Since packing in bytes, byteorder is mostly irrelevant
            # except when verifying zero padding at end of file
            (data,) = unpack("<I", word)  # comma for unpacking tuple result
            data_b = data.to_bytes(4, byteorder="big")
            detections.extend(data_b)

    # Verify detections match length entry in header, ignoring \x00 byte padding
    if headt3.length_entry != len(detections) - data_b.count(0):
        logger.error("Stream 3 size inconsistency: length of entry does not match.")

    # Assign to coincidence matrix
    for i in range(headt3.length_entry):
        detection = detections[i]
        a = decode[(detection >> 4) & 0xF]  # Alice
        b = decode[detection & 0xF]  # Bob
        if a >= 0 and b >= 0:
            service.coinc_matrix[a * 4 + b] += 1
            service.okcount += 1
        else:
            if a < 0:
                service.garbage[0] += 1
            if b < 0:
                service.garbage[1] += 1

    # Calculate QBER
    er_coin = sum(service.coinc_matrix[i] for i in er_coinc_id)
    gd_coin = sum(service.coinc_matrix[i] for i in gd_coinc_id)
    if er_coin + gd_coin != 0:
        service.qber = float(round(er_coin / (er_coin + gd_coin), 3))  # ignore garbage

    return service


def read_T4(filename: str):
    """Returns timestamp indices and detector information from T4 files.

    Timestamp and detector information follow the formats specified in
    'parse_timestamps' module, for interoperability.

    Raises:
        ValueError: Length and termination bits inconsistent with filespec.
    """

    # Header details
    header = read_T4_header(filename)
    timeorder = header.timeorder
    timeorder_extended = 32
    basebits = header.base_bits
    length = header.length

    # Accumulators
    indices = []
    bases = []
    time_index = 0
    buffer = 0
    size = 0  # number of bits in buffer
    with open(filename, "rb") as f:
        # f = BufferedFileRead(f)
        f.read(20)  # remove header
        while True:
            # Read timing information
            time_dindex, buffer, size = extract_bits(timeorder, buffer, size, f)

            # End-of-file check
            if time_dindex == 1:
                base, buffer, size = extract_bits(basebits, buffer, size, f)
                if base != 0:
                    raise ValueError(
                        "File inconsistent with T4 epoch definition: "
                        "Terminal base bits not zero"
                    )
                break

            # Check if timing is actually in extended format, i.e. 32-bits
            if time_dindex == 0:
                time_dindex, buffer, size = extract_bits(
                    timeorder_extended, buffer, size, f
                )

            # Index retrieved, normalize
            time_index += time_dindex - 2
            indices.append(time_index)

            # Extract detector pattern, but not important here
            # Note we are guaranteed (timeorder+basebits < 32)
            base, buffer, size = extract_bits(basebits, buffer, size, f)
            bases.append(base)

    # Validity checks
    if length and len(indices) != length:
        raise ValueError("Number of epochs do not match length specified in header")

    return np.array(indices), np.array(bases)
