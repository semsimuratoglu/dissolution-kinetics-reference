"""Goodness-of-fit metrics computed in original percentage-dissolved space."""

from __future__ import annotations

import numpy as np


def rss(observed: np.ndarray, predicted: np.ndarray) -> float:
    """Residual sum of squares."""
    residual = observed - predicted
    return float(np.sum(residual ** 2))


def tss(observed: np.ndarray) -> float:
    """Total sum of squares about the mean of observed values."""
    mean_obs = float(np.mean(observed))
    return float(np.sum((observed - mean_obs) ** 2))


def r_squared(observed: np.ndarray, predicted: np.ndarray) -> float:
    """R² = 1 - RSS / TSS."""
    total = tss(observed)
    if total == 0:
        return 1.0 if rss(observed, predicted) == 0 else 0.0
    return 1.0 - rss(observed, predicted) / total


def adjusted_r_squared(
    observed: np.ndarray, predicted: np.ndarray, n_params: int
) -> float | None:
    """adj-R² = 1 - (1 - R²) * (N - 1) / (N - P - 1).

    Returns ``None`` when ``N - P - 1 <= 0`` (insufficient degrees of freedom).
    """
    n = len(observed)
    denom = n - n_params - 1
    if denom <= 0:
        return None
    r2 = r_squared(observed, predicted)
    return 1.0 - (1.0 - r2) * (n - 1) / denom


def rmse(observed: np.ndarray, predicted: np.ndarray, n_params: int) -> float | None:
    """RMSE = sqrt(RSS / (N - P))  (DDSolver convention, DF = n - p)."""
    n = len(observed)
    df = n - n_params
    if df <= 0:
        return None
    return float(np.sqrt(rss(observed, predicted) / df))


def aic(observed: np.ndarray, predicted: np.ndarray, n_params: int) -> float:
    """AIC = N * ln(RSS) + 2 * P  (DDSolver convention).

    Uses the natural logarithm. RSS must be strictly positive; a small floor is
    applied only when RSS is numerically zero to avoid ``log(0)``.
    """
    n = len(observed)
    residual_sum = rss(observed, predicted)
    if residual_sum <= 0:
        residual_sum = np.finfo(float).tiny
    return float(n * np.log(residual_sum) + 2 * n_params)


def msc(observed: np.ndarray, predicted: np.ndarray, n_params: int) -> float:
    """MSC = ln(TSS / RSS) - 2 * P / N  (model selection criterion)."""
    n = len(observed)
    total = tss(observed)
    residual_sum = rss(observed, predicted)
    if residual_sum <= 0:
        residual_sum = np.finfo(float).tiny
    if total <= 0:
        total = np.finfo(float).tiny
    return float(np.log(total / residual_sum) - 2 * n_params / n)


def goodness_of_fit(
    observed: np.ndarray, predicted: np.ndarray, n_params: int
) -> dict[str, float | None]:
    """Compute all GOF metrics for one fit."""
    return {
        "R2": r_squared(observed, predicted),
        "adj_R2": adjusted_r_squared(observed, predicted, n_params),
        "RMSE": rmse(observed, predicted, n_params),
        "AIC": aic(observed, predicted, n_params),
        "MSC": msc(observed, predicted, n_params),
        "RSS": rss(observed, predicted),
    }
