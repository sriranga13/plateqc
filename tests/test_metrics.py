import math

import pytest

from plateqc.metrics import (
    mean, stdev, cv, zprime, ssmd, signal_window, s2b,
    edge_ratio, drift, row_cv, col_cv,
)

HIGH = [1000.0] * 8
LOW = [100.0] * 8


def test_mean_empty():
    assert mean([]) is None


def test_stdev_needs_two():
    assert stdev([5.0]) is None
    assert stdev([]) is None
    assert stdev([1.0, 3.0]) == pytest.approx(math.sqrt(2.0))


def test_cv():
    assert cv([100.0, 100.0, 100.0]) == pytest.approx(0.0)
    assert cv([]) is None
    assert cv([0.0, 0.0]) is None


def test_zprime_perfect_separation():
    assert zprime(HIGH, LOW) == pytest.approx(1.0)


def test_zprime_typical():
    high = [1000, 980, 1020, 990, 1010, 995, 1005, 985]
    low = [100, 95, 105, 98, 102, 99, 101, 103]
    z = zprime(high, low)
    assert 0.5 < z < 1.0


def test_zprime_no_separation():
    assert zprime(HIGH, HIGH) is None  # zero denominator


def test_zprime_needs_variance():
    assert zprime([1.0], LOW) is None


def test_ssmd():
    high = [1000, 980, 1020, 990, 1010, 995, 1005, 985]
    low = [100, 95, 105, 98, 102, 99, 101, 103]
    s = ssmd(high, low)
    assert s > 5.0


def test_signal_window():
    assert signal_window(HIGH, LOW) == pytest.approx(900.0)
    assert signal_window([], LOW) is None


def test_s2b():
    assert s2b(HIGH, LOW) == pytest.approx(10.0)
    assert s2b(HIGH, [0.0] * 8) is None


def test_edge_ratio():
    assert edge_ratio([110.0] * 10, [100.0] * 50) == pytest.approx(0.1)
    assert edge_ratio([], [100.0]) is None


def test_drift_flat():
    assert drift([0, 1, 2, 3], [100.0] * 4) == pytest.approx(0.0)


def test_drift_ten_percent():
    # rises from 100 to 110 across 4 columns -> ~10% drift
    assert drift([0, 1, 2, 3], [100.0, 103.33, 106.67, 110.0]) == pytest.approx(0.1, abs=1e-2)


def test_drift_constant_x():
    assert drift([1, 1, 1], [1.0, 2.0, 3.0]) is None


def test_row_col_cv():
    assert row_cv([100.0, 100.0, 100.0]) == pytest.approx(0.0)
    assert col_cv([]) is None
