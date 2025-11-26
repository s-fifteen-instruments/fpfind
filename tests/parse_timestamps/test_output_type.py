import numpy as np
import pytest

from fpfind import NP_PRECISEFLOAT, TSRES
from fpfind.lib import parse_timestamps as parser

testcases = [
    (TSRES.NS1, False, np.uint64),
    (TSRES.NS1, True, NP_PRECISEFLOAT),
    (TSRES.NS2, False, np.uint64),
    (TSRES.NS2, True, NP_PRECISEFLOAT),
    (TSRES.PS125, False, np.uint64),
    (TSRES.PS125, True, NP_PRECISEFLOAT),
    (TSRES.PS4, False, np.uint64),
    (TSRES.PS4, True, NP_PRECISEFLOAT),
]


@pytest.mark.parametrize("resolution,fractional,rtype", testcases)
def test_reada0_type(path_timestamp_a0, resolution, fractional, rtype):
    t, p = parser.read_a0(
        path_timestamp_a0,
        resolution=resolution,
        fractional=fractional,
    )
    assert t.dtype == rtype
    assert p.dtype == np.uint32


@pytest.mark.parametrize("resolution,fractional,rtype", testcases)
def test_reada1_type(path_timestamp_a1, resolution, fractional, rtype):
    t, p = parser.read_a1(
        path_timestamp_a1,
        legacy=True,
        resolution=resolution,
        fractional=fractional,
    )
    assert t.dtype == rtype
    assert p.dtype == np.uint32


@pytest.mark.parametrize("resolution,fractional,rtype", testcases)
def test_reada2_type(path_timestamp_a2, resolution, fractional, rtype):
    t, p = parser.read_a2(
        path_timestamp_a2,
        resolution=resolution,
        fractional=fractional,
    )
    assert t.dtype == rtype
    assert p.dtype == np.uint32
