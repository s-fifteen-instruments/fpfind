import numpy as np

from fpfind import TSRES
from fpfind.lib import parse_timestamps as parser


def test_reada1_1nsint(path_timestamp_a1):
    t, p = parser.read_a1(
        path_timestamp_a1,
        legacy=True,
        resolution=TSRES.NS1,
        fractional=False,
    )
    assert t.dtype == np.uint64
    assert p.dtype == np.uint32


def test_reada1_1nsfloat(path_timestamp_a1):
    t, p = parser.read_a1(
        path_timestamp_a1,
        legacy=True,
        resolution=TSRES.NS1,
        fractional=True,
    )
    assert t.dtype == np.float128
    assert p.dtype == np.uint32


def test_reada1_4psint(path_timestamp_a1):
    t, p = parser.read_a1(
        path_timestamp_a1,
        legacy=True,
        resolution=TSRES.PS4,
        fractional=False,
    )
    assert t.dtype == np.uint64
    assert p.dtype == np.uint32


def test_reada1_4psfloat(path_timestamp_a1):
    t, p = parser.read_a1(
        path_timestamp_a1,
        legacy=True,
        resolution=TSRES.PS4,
        fractional=True,
    )
    assert t.dtype == np.float128
    assert p.dtype == np.uint32
