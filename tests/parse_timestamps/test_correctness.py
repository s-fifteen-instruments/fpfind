import numpy as np

import fpfind.lib.parse_timestamps as parser
from fpfind import TSRES


def test_read_a0(path_timestamp_a0):
    ts, ps = parser.read_a0(path_timestamp_a0)
    assert ts.size == ps.size == 1000
    assert np.unique(ps).tolist() == [1, 8]


def test_read_a1(path_timestamp_a1):
    ts, ps = parser.read_a1(path_timestamp_a1, legacy=True)
    assert ts.size == ps.size == 1000
    assert np.unique(ps).tolist() == [1, 8]


def test_read_a1_rollover(path_timestamp_a1):
    ts, ps = parser.read_a1(path_timestamp_a1, legacy=True, ignore_rollover=True)
    assert ts.size == ps.size == 862
    assert np.unique(ps).tolist() == [1, 8]

    # Check equal number of events
    ts1 = ts[(ps & 0x1).astype(bool)]
    ts2 = ts[(ps & 0x8).astype(bool)]
    assert ts1.size == ts2.size


def test_read_equivalence_1ns(path_timestamp_a0, path_timestamp_a1, path_timestamp_a2):
    # Equivalence up to 2**-54 = 5.55e-17
    rtol = 6e-17

    ts0, ps0 = parser.read_a0(path_timestamp_a0)
    ts1, ps1 = parser.read_a1(path_timestamp_a1, legacy=True)
    ts2, ps2 = parser.read_a2(path_timestamp_a2)
    assert np.allclose(ts0, ts1, rtol=rtol, atol=0)
    assert np.allclose(ts0, ts2, rtol=rtol, atol=0)
    assert np.all(ps0 == ps1)
    assert np.all(ps0 == ps2)


def test_read_equivalence_4ps(path_timestamp_a0, path_timestamp_a1, path_timestamp_a2):
    # Verifies timestamp and detector pattern values are exactly the same
    kwargs = {
        "resolution": TSRES.PS4,
        "fractional": False,
    }
    ts0, ps0 = parser.read_a0(path_timestamp_a0, **kwargs)
    ts1, ps1 = parser.read_a1(path_timestamp_a1, **kwargs, legacy=True)
    ts2, ps2 = parser.read_a2(path_timestamp_a2, **kwargs)
    assert np.all(ts0 == ts1)
    assert np.all(ts0 == ts2)
    assert np.all(ps0 == ps1)
    assert np.all(ps0 == ps2)


def test_read_equivalence_4ps_rollover(
    path_timestamp_a0, path_timestamp_a1, path_timestamp_a2
):
    # Verifies timestamp and detector pattern values are exactly the same,
    # even after filtering out dummy rollover events
    kwargs = {
        "resolution": TSRES.PS4,
        "fractional": False,
        "ignore_rollover": True,
    }
    ts0, ps0 = parser.read_a0(path_timestamp_a0, **kwargs)
    ts1, ps1 = parser.read_a1(path_timestamp_a1, **kwargs, legacy=True)
    ts2, ps2 = parser.read_a2(path_timestamp_a2, **kwargs)
    assert np.all(ts0 == ts1)
    assert np.all(ts0 == ts2)
    assert np.all(ps0 == ps1)
    assert np.all(ps0 == ps2)
