# Test different datasets with fine-tuned configuration

import subprocess

import pytest

datadir = "tests/data/internal"

fast_commands = [
    f"-t {datadir}/1511_qkd/20251107_dataset304a1_noref_noref_alice_vis.ts "
    f"-T {datadir}/1511_qkd/20251107_dataset304a1_noref_noref_alice_vis.ts "
    "-q19 -R16 --df=-9e-6",
]

slow_commands = [

]

@pytest.mark.parametrize("command", fast_commands)
@pytest.mark.internal
def test_fpfind_functional_fast(command):
    subprocess.check_call(f"fpfind {command}", shell=True)


@pytest.mark.parametrize("command", slow_commands)
@pytest.mark.internal
@pytest.mark.slow
def test_fpfind_functional_slow(command):
    pass
