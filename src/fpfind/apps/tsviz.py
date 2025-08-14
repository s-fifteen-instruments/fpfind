#!/usr/bin/env python3
# Visualize timestamp file
# Justin, 2025-04-07
#
# Currently hardcoded with read_a1(legacy=True).

import sys

import matplotlib.pyplot as plt
import numpy as np

import fpfind.lib.parse_timestamps as parser

if len(sys.argv) < 2:
    print("usage: tsviz.py <FILE> [legacy]")
    sys.exit(1)

legacy = len(sys.argv) >= 3
filename = sys.argv[1]
ts, ps = parser.read_a1(filename, legacy=legacy, ignore_rollover=True)


def generate_xy(ts, ps, ch):
    """
    Args:
        ch: Channel number, i.e. between 1 and 4 inclusive.
    """
    p = int(2 ** (ch - 1))
    xs = ts[(ps & p).astype(bool)]
    ys = np.ones(xs.size) * ch
    return xs, ys


args = ["x"]
kwargs = {"markersize": 3}

fig, ax = plt.subplots(figsize=(12, 3))
plt.plot(*generate_xy(ts, ps, 1), *args, **kwargs)
plt.plot(*generate_xy(ts, ps, 2), *args, **kwargs)
plt.plot(*generate_xy(ts, ps, 3), *args, **kwargs)
plt.plot(*generate_xy(ts, ps, 4), *args, **kwargs)

plt.ylim([-1, 6])
plt.yticks([1, 2, 3, 4])
plt.gca().invert_yaxis()
plt.ylabel("Channel")
plt.xlabel("Time (ns)")
plt.tight_layout()

plt.show()
