from __future__ import annotations

import numpy as np
import pytest

from panacea.utils import safe_division, build_weight_matrix, robust_polyfit
from panacea.routine import fit_response_cont


def test_safe_division_basic_and_masking():
    # Same-shape arrays: zeros/inf in denom should yield fill value
    num = np.array([1.0, 2.0, -3.0, 4.0])
    denom = np.array([1.0, 0.0, np.inf, 1e-12])
    out = safe_division(num, denom, eps=1e-8, fillval=-999.0)
    assert np.isclose(out[0], 1.0)
    # denom == 0 -> fill
    assert np.isclose(out[1], -999.0)
    # denom == inf -> 0 (treated as good but  num/inf -> 0); however function flags finite only, so inf considered bad
    assert np.isclose(out[2], -999.0)
    # tiny denom <= eps -> fill
    assert np.isclose(out[3], -999.0)


def test_safe_division_broadcast_over_columns():
    # 2D numerator with 1D denom triggers broadcasting branch
    rng = np.random.default_rng(123)
    num = rng.normal(size=(4, 6))
    denom = np.array([1.0, 2.0, 0.0, 4.0, 1e-12, np.nan])
    out = safe_division(num, denom, eps=1e-8, fillval=0.0)
    # Good columns 0,1,3 should equal naive division
    for j in [0, 1, 3]:
        assert np.allclose(out[:, j], num[:, j] / denom[j])
    # Bad columns 2,4,5 should be fillval
    for j in [2, 4, 5]:
        assert np.all(out[:, j] == 0.0)


def test_build_weight_matrix_properties():
    # Arrange fibers on a small grid
    x = np.array([0.0, 1.0, 0.0, 1.0])
    y = np.array([0.0, 0.0, 1.0, 1.0])
    W = build_weight_matrix(x, y, sig=0.75)
    # Column normalization: columns sum to 1.0
    colsum = W.sum(axis=0)
    assert np.allclose(colsum, 1.0, atol=1e-6)
    # Self-weight should be the largest in each column
    diag = np.diag(W)
    assert np.all(diag >= (W.max(axis=0) - 1e-12))
    # Non-negativity
    assert np.all(W >= 0)


@pytest.mark.filterwarnings("ignore:polyfit may be poorly conditioned")
def test_robust_polyfit_handles_outliers_and_noise():
    rng = np.random.default_rng(42)
    x = np.linspace(-1, 1, 201)
    # True quadratic
    y_true = 0.5 * x**2 - 0.2 * x + 1.0
    y = y_true + rng.normal(scale=0.03, size=x.size)
    # Inject some large outliers
    idx = rng.choice(x.size, size=10, replace=False)
    y[idx] += rng.normal(loc=0.0, scale=1.5, size=idx.size)
    # Ensure positives to satisfy initial mask y>0.0
    y = np.abs(y)
    y_fit = robust_polyfit(x, y, order=2, niter=5)
    # Evaluate RMS error; should be small compared to outlier scale
    rms = np.sqrt(np.mean((y_fit - y_true) ** 2))
    assert rms < 0.08


def test_fit_response_cont_smooths_and_preserves_shape():
    # Construct a synthetic sky spectrum: smooth continuum + sharp lines
    rng = np.random.default_rng(7)
    wv = np.linspace(4000.0, 5000.0, 501)
    cont_true = 1.0 + 0.001 * (wv - wv.mean())  # gentle slope
    sky = cont_true.copy()
    # Add narrow emission/absorption features
    line_centers = [4050, 4200, 4500, 4700, 4950]
    for c in line_centers:
        j = np.argmin(np.abs(wv - c))
        sky[j] += 0.8  # emission spike
        if j + 2 < sky.size:
            sky[j + 2] -= 0.5  # nearby absorption
    # Add small noise
    sky += rng.normal(scale=0.02, size=sky.size)

    cont_est = fit_response_cont(wv, sky, skip=4, fil_len=21, func=np.array)

    # Sanity checks
    assert cont_est.shape == sky.shape
    # The estimate should be close to the true continuum away from edges
    core = slice(20, -20)
    mae = np.mean(np.abs(cont_est[core] - cont_true[core]))
    assert mae < 0.05
    # Smoothing should reduce variance relative to original
    assert np.var(cont_est[core] - np.mean(cont_est[core])) < np.var(sky[core] - np.mean(sky[core]))
