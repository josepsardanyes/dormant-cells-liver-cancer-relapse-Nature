"""Estimate the untreated exponential tumour-growth rate."""

import numpy as np
from scipy.stats import linregress


def fit_gs(time, volume):
    """Return growth-rate slope and R-squared for log(volume) versus time."""
    time = np.asarray(time, dtype=float)
    volume = np.asarray(volume, dtype=float)
    if time.shape != volume.shape:
        raise ValueError("time and volume must have the same shape")

    valid = np.isfinite(time) & np.isfinite(volume) & (volume > 0)
    if np.count_nonzero(valid) < 2:
        raise ValueError("at least two finite observations with positive volume are required")

    regression = linregress(time[valid], np.log(volume[valid]))
    return regression.slope, regression.rvalue**2
