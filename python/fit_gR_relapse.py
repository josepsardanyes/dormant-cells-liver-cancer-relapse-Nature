"""Estimate the exponential tumour-growth rate during relapse."""

import numpy as np
from scipy.stats import linregress


def fit_gr(time, volume, start_day=26):
    """Fit log(volume) from the post-start-day volume minimum onward."""
    time = np.asarray(time, dtype=float)
    volume = np.asarray(volume, dtype=float)
    if time.shape != volume.shape:
        raise ValueError("time and volume must have the same shape")

    valid = np.isfinite(time) & np.isfinite(volume) & (volume > 0)
    t = time[valid]
    v = volume[valid]
    post = np.flatnonzero(t >= start_day)
    if post.size == 0:
        raise ValueError("no valid observations occur on or after start_day")

    index_minimum = post[np.argmin(v[post])]
    if t[index_minimum:].size < 2:
        raise ValueError("at least two observations are required from the relapse minimum onward")

    regression = linregress(t[index_minimum:], np.log(v[index_minimum:]))
    return regression.slope, regression.rvalue**2
