#!/usr/bin/env python3
# Check correctness of timestamp outputs

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit

from fpfind.lib.parse_timestamps import read_a1

reference = "data/sample_ts.dat"
target = "data/sample_ts.4000ppb.dat"
drift = 3999.9722503125668

# reference = "test_10s2"
# target = "test_10s2_1000ppb"
# drift = 1000

ts, ps = read_a1(reference, legacy=True)
ts2, ps2 = read_a1(target, legacy=True)

error = ts2 - ts
ts0 = (ts - ts[0]) / 1e9  # shifted to zero in [s]
expected = ts0 * drift


# Estimate the slope
def linear(t, m):
    return m * t


popt, pcov = curve_fit(linear, ts0, error, p0=(drift,))
perr = np.sqrt(np.diag(pcov))
print(popt, perr)
print(error[-1] - expected[-1])
plt.plot(ts0, error, linewidth=1, label="error")
plt.plot(ts0, expected, "x", markersize=3, label="expected")
plt.plot(ts0, linear(ts0, *popt), "--", linewidth=1, label="estimated")
# plt.xlim(left=15.931, right=15.932)
# plt.ylim([63720,63730])
plt.legend()
plt.show()

"""
>>> 2**-34 * 1e9
0.05820766091346741  # resolution in ppb

>>> 4000 / (2**-34 * 1e9)
68719.476736  # in units of 2**-34

>>> 68719 * 2**-34 * 1e9
3999.9722503125668  # expected drift
"""
