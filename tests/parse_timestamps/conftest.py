import pytest


@pytest.fixture
def path_timestamp_a0(path_tsdata):
    return path_tsdata / "timestamp7_c14_rollover.Aa0f.ts"


@pytest.fixture
def path_timestamp_a1(path_tsdata):
    return path_tsdata / "timestamp7_c14_rollover.Aa1fX.ts"


@pytest.fixture
def path_timestamp_a2(path_tsdata):
    return path_tsdata / "timestamp7_c14_rollover.Aa2f.ts"
