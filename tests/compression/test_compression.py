import shutil

import numpy as np
import pytest

from fpfind.lib.compression import DTS
import fpfind.lib.parse_timestamps as parser

LOAD_RAW = {
    "resolution": parser.TSRES.PS4,
    "fractional": False,
}

@pytest.fixture
def path_timestamp(path_tsdata):
    return path_tsdata / "timestamp7_c1234.Aa1.ts"

def test_DTS_load_decoded_lossless_internal(path_timestamp):
    ts, ps = parser.read_a1(path_timestamp, **LOAD_RAW)
    dts = DTS._load_decoded(path_timestamp, lossless=True)
    _ts, _ps = dts._get()
    _ps &= 0xF  # extract only detector information
    assert np.all(ts == _ts)
    assert np.all(ps == _ps)

def test_DTS_load_decoded_internal(path_timestamp):
    ts, ps = parser.read_a1(path_timestamp, **LOAD_RAW)
    dts = DTS._load_decoded(path_timestamp)
    _ts, _ps = dts._get()
    assert np.all(ts == _ts)
    assert np.all(ps == _ps)

def test_DTS_decoded_roundtrip_lossless_internal(path_timestamp, tmp_path):
    original = tmp_path / "original.ts"
    restored = tmp_path / "restored.ts"
    shutil.copy(path_timestamp, original, follow_symlinks=True)

    dts = DTS._load_decoded(original, lossless=True)
    dts._save_decoded(restored)

    ts, ps = parser.read_a1(original, **LOAD_RAW)
    _ts, _ps = parser.read_a1(restored, **LOAD_RAW)
    assert np.all(ts == _ts)
    assert np.all(ps == _ps)

def test_DTS_decoded_roundtrip_internal(path_timestamp, tmp_path):
    original = tmp_path / "original.ts"
    restored = tmp_path / "restored.ts"
    shutil.copy(path_timestamp, original, follow_symlinks=True)

    dts = DTS._load_decoded(original)
    dts._save_decoded(restored)

    ts, ps = parser.read_a1(original, resolution=parser.TSRES.PS4, fractional=False)
    _ts, _ps = parser.read_a1(restored, resolution=parser.TSRES.PS4, fractional=False)
    assert np.all(ts == _ts)
    assert np.all(ps == _ps)

def test_DTS_roundtrip_lossless(path_timestamp, tmp_path):
    original = tmp_path / "original.ts"
    encoded = original.with_suffix(".dts")
    restored = tmp_path / "restored.ts"
    shutil.copy(path_timestamp, original, follow_symlinks=True)

    DTS.encode(original, lossless=True, target=encoded)
    DTS.decode(encoded, target=restored)

    ts, ps = parser.read_a1(original, **LOAD_RAW)
    _ts, _ps = parser.read_a1(restored, **LOAD_RAW)
    assert np.all(ts == _ts)
    assert np.all(ps == _ps)

    # Test compression performance
    assert encoded.stat().st_size < original.stat().st_size * 0.70

def test_DTS_roundtrip(path_timestamp, tmp_path):
    original = tmp_path / "original.ts"
    encoded = original.with_suffix(".dts")
    restored = tmp_path / "restored.ts"
    shutil.copy(path_timestamp, original, follow_symlinks=True)

    DTS.encode(original, target=encoded)
    DTS.decode(encoded, target=restored)

    ts, ps = parser.read_a1(original, **LOAD_RAW)
    _ts, _ps = parser.read_a1(restored, **LOAD_RAW)
    assert np.all(ts == _ts)
    assert np.all(ps == _ps)

    # Test compression performance
    assert encoded.stat().st_size < original.stat().st_size * 0.55

def test_DTS_load_lossless(path_timestamp, tmp_path):
    original = tmp_path / "original.ts"
    encoded = original.with_suffix(".dts")
    shutil.copy(path_timestamp, original, follow_symlinks=True)

    DTS.encode(original, lossless=True, target=encoded)

    ts, ps = parser.read_a1(original, **LOAD_RAW)
    _ts, _ps = DTS.load(encoded, **LOAD_RAW)
    assert np.all(ts == _ts)
    assert np.all(ps == _ps)

def test_DTS_load(path_timestamp, tmp_path):
    original = tmp_path / "original.ts"
    encoded = original.with_suffix(".dts")
    shutil.copy(path_timestamp, original, follow_symlinks=True)

    DTS.encode(original, target=encoded)

    ts, ps = parser.read_a1(original, **LOAD_RAW)
    _ts, _ps = DTS.load(encoded, **LOAD_RAW)
    assert np.all(ts == _ts)
    assert np.all(ps == _ps)
