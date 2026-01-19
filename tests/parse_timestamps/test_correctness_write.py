import struct

import numpy as np

import fpfind.lib.parse_timestamps as parser
from fpfind import TSRES


def _write_a1_deprec(
    filename,
    t,
    p,
    legacy: bool = False,
    resolution: TSRES = TSRES.NS1,
    allow_duplicates: bool = True,
) -> None:
    """Ported from 'fpfind.lib.parse_timestamps.write_a1()'.

    Preserved here for correctness testing.
    """
    events = parser._consolidate_events(t, p, resolution, allow_duplicates=allow_duplicates)
    with open(filename, "wb", opener=parser.stdopen) as f:
        for line in events:
            if legacy:
                line = int(line)
                line = ((line & 0xFFFFFFFF) << 32) + (line >> 32)
            f.write(struct.pack("=Q", line))


def test_write_a1(path_timestamp_a1, tmp_path):
    target = tmp_path / "ts"
    ts, ps = parser.read_a1(path_timestamp_a1, legacy=True)

    parser.write_a1(target, ts, ps)
    _ts, _ps = parser.read_a1(target)
    assert np.all(ts == _ts)
    assert np.all(ps == _ps)

    # With legacy mode
    parser.write_a1(target, ts, ps, legacy=True)
    _ts, _ps = parser.read_a1(target, legacy=True)
    assert np.all(ts == _ts)
    assert np.all(ps == _ps)


def test_write_a1_deprec(path_timestamp_a1, tmp_path):
    target = tmp_path / "ts"
    target_deprec = tmp_path / "ts.deprec"
    ts, ps = parser.read_a1(path_timestamp_a1, legacy=True)

    parser.write_a1(target, ts, ps)
    _write_a1_deprec(target_deprec, ts, ps)
    ts, ps = parser.read_a1(target)
    _ts, _ps = parser.read_a1(target_deprec)
    assert np.all(ts == _ts)
    assert np.all(ps == _ps)
