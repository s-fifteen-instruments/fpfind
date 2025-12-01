# Test with same dataset across different configuration

import subprocess

import pytest


@pytest.fixture
def alice(intdatadir):
    return f"{intdatadir}/1511_qkd/20251107_dataset304a1_noref_noref_alice_vis.ts"

@pytest.fixture
def bob(intdatadir):
    return f"{intdatadir}/1511_qkd/20251107_dataset304a2_noref_noref_bob_ir.ts"


options = [
    "-q19 -R16 --df=-9e-6",
    "-q19 -R16 --df=-8e-6",
    "-q19 -R16 --df=-7e-6",
    "-q20 -R8 --df=-7e-6",
]

@pytest.mark.parametrize("opts", options)
@pytest.mark.internal
def test_fpfind(alice, bob, opts):
    subprocess.check_call(f"fpfind -t {alice} -T {bob} {opts}", shell=True)
