import sys
import warnings

import numpy as np
import pytest

# Automatically skip 80-bit precision checks if not supported
# Note: np.double only goes up to 1.79e308, while np.longdouble
#       supports up to 1.19e4932 (if implemented as 80-bit float).
#       For more information, see 'np.finfo(np.longdouble)'.
is_extfloat_supported = np.longdouble("1e309") != np.inf
if not is_extfloat_supported:
    warnings.warn("Extended precision float (80-bit) is not supported.")


def pytest_addoption(parser):
    parser.addoption(
        "--internal", action="store_true", default=False, help="run internal tests"
    )
    parser.addoption(
        "--slow", action="store_true", default=False, help="run slow tests"
    )


def pytest_collection_modifyitems(config, items):
    def skip(marker, reason):
        for item in items:
            if marker in item.keywords:
                skip_marker = pytest.mark.skip(reason=reason)
                item.add_marker(skip_marker)

    if not config.getoption("--internal"):
        skip("internal", "need --internal")

    if not config.getoption("--slow"):
        skip("slow", "need --slow")

    if is_extfloat_supported:
        skip("float64", "extended precision supported")
    else:
        skip("float80", "no support for extended precision")

    # Exclude Windows from Linux-only tests
    if sys.platform == "win32":
        skip("unix", "cannot run on Windows")


@pytest.fixture
def testdir(pytestconfig):
    return pytestconfig.rootpath / "tests"


@pytest.fixture
def datadir(testdir):
    return testdir / "data"


@pytest.fixture
def intdatadir(datadir):
    return datadir / "internal"


# To deprecate
@pytest.fixture
def path_tests(testdir):
    return testdir


# To deprecate
@pytest.fixture
def path_tsdata(datadir):
    return datadir
