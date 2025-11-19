from importlib.metadata import version

import fpfind.lib
from fpfind.lib.constants import NP_PRECISEFLOAT, TSRES
from fpfind.lib.utils import get_overlap

__all__ = [
    "fpfind",
    "NP_PRECISEFLOAT",
    "TSRES",
    "get_overlap",  # used in tests
]

VERSION = version("fpfind")
