"""Assay QC statistics: Z'-factor, SSMD, signal window, edge/drift diagnostics.

All metrics follow Zhang et al. (J Biomol Screen 1999) conventions:
Z' = 1 - 3*(sd_high + sd_low) / |mean_high - mean_low|.
Sample (n-1) standard deviation is used; metrics return None when the
inputs cannot support them instead of raising.
"""

import math


def mean(xs):
    xs = list(xs)
    if not xs:
        return None
    return sum(xs) / len(xs)


def stdev(xs):
    """Sample standard deviation; None for fewer than 2 values."""
    xs = list(xs)
    n = len(xs)
    if n < 2:
        return None
    m = mean(xs)
    return math.sqrt(sum((x - m) ** 2 for x in xs) / (n - 1))


def cv(xs):
    """Coefficient of variation (sd/mean); None if undefined."""
    m, s = mean(xs), stdev(xs)
    if m is None or s is None or m == 0:
        return None
    return s / abs(m)


def zprime(high, low):
    """Z'-factor for high vs low control wells. 1 is ideal, >0.5 excellent."""
    mh, ml = mean(high), mean(low)
    sh, sl = stdev(high), stdev(low)
    if None in (mh, ml, sh, sl):
        return None
    denom = abs(mh - ml)
    if denom == 0:
        return None
    return 1 - (3 * (sh + sl)) / denom


def ssmd(high, low):
    """Strictly standardized mean difference, high vs low controls."""
    mh, ml = mean(high), mean(low)
    sh, sl = stdev(high), stdev(low)
    if None in (mh, ml, sh, sl):
        return None
    pooled = math.sqrt(sh ** 2 + sl ** 2)
    if pooled == 0:
        return None
    return abs(mh - ml) / pooled


def signal_window(high, low):
    """Separation band: (mh - 3*sh) - (ml + 3*sl); positive is good."""
    mh, ml = mean(high), mean(low)
    sh, sl = stdev(high), stdev(low)
    if None in (mh, ml, sh, sl):
        return None
    return (mh - 3 * sh) - (ml + 3 * sl)


def s2b(high, low):
    """Signal-to-background ratio, mean(high)/mean(low)."""
    mh, ml = mean(high), mean(low)
    if mh is None or ml is None or ml == 0:
        return None
    return mh / ml


def edge_ratio(edge_values, interior_values):
    """mean(edge)/mean(interior) - 1; positive means edge wells read high.

    A classic evaporation/temperature edge effect shows up here before it
    is visible in a heatmap.
    """
    me, mi = mean(edge_values), mean(interior_values)
    if me is None or mi is None or mi == 0:
        return None
    return me / mi - 1


def _linregress(xs, ys):
    """Return (slope, intercept); None, None when undefined."""
    xs, ys = list(xs), list(ys)
    n = len(xs)
    if n < 2 or len(set(xs)) < 2:
        return None, None
    mx, my = mean(xs), mean(ys)
    denom = sum((x - mx) ** 2 for x in xs)
    if denom == 0:
        return None, None
    slope = sum((x - mx) * (y - my) for x, y in zip(xs, ys)) / denom
    return slope, my - slope * mx


def drift(column_indices, values):
    """Fractional drift per column: slope * span / mean.

    E.g. 0.10 means the signal drifts ~10% across the plate's columns,
    hinting at dispense timing or reader warm-up effects.
    """
    xs, ys = list(column_indices), list(values)
    m = mean(ys)
    slope, _ = _linregress(xs, ys)
    if m is None or slope is None or m == 0 or len(set(xs)) < 2:
        return None
    span = max(xs) - min(xs)
    return slope * span / abs(m)


def row_cv(row_means):
    """CV of per-row means: flags rows that behave differently."""
    return cv(row_means)


def col_cv(col_means):
    """CV of per-column means: flags columns that behave differently."""
    return cv(col_means)
