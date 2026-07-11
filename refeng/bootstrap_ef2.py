"""Extended E(f2) bootstrap: percentile + BCa CIs, granularity comparison, validation mode."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.stats import norm

from refeng.metrics import expected_f2_analytical, f2, fda_truncation, ema_truncation


@dataclass
class BootstrapEf2Validation:
    """Full bootstrap validation output for external anchoring."""

    f2_observed: float
    e_f2_point: float
    bc_f2: float
    percentile_ci_lower: float
    percentile_ci_upper: float
    bca_ci_lower: float
    bca_ci_upper: float
    median: float
    n_iterations: int
    seed: int
    truncation: str
    similar_percentile: bool
    similar_bca: bool


def _primary_truncation_rule(times: np.ndarray, ref_means: np.ndarray, test_means: np.ndarray) -> str:
    fda = fda_truncation(times, ref_means, test_means)
    ema = ema_truncation(times, ref_means, test_means)
    if np.array_equal(fda.kept_indices, ema.kept_indices):
        return "FDA"
    return "EMA"


def _observed_f2_primary(times: np.ndarray, ref_means: np.ndarray, test_means: np.ndarray) -> float:
    rule = _primary_truncation_rule(times, ref_means, test_means)
    trunc = ema_truncation(times, ref_means, test_means) if rule == "EMA" else fda_truncation(times, ref_means, test_means)
    return f2(trunc.ref_means, trunc.test_means)


def _resolve_trunc_rule(times: np.ndarray, ref_units: np.ndarray, test_units: np.ndarray, truncation: str) -> str:
    ref_means = np.mean(ref_units, axis=0)
    test_means = np.mean(test_units, axis=0)
    rule = _primary_truncation_rule(times, ref_means, test_means)
    if truncation.upper() == "EMA" or rule == "EMA":
        return "EMA"
    return "FDA"


def _jackknife_ef2_values(
    times: np.ndarray,
    ref_units: np.ndarray,
    test_units: np.ndarray,
    trunc_rule: str,
) -> np.ndarray:
    """Delete-one jackknife estimates of E(f2) for BCa acceleration."""
    n_ref, n_test = ref_units.shape[0], test_units.shape[0]
    vals: list[float] = []
    for i in range(n_ref):
        v = expected_f2_analytical(times, np.delete(ref_units, i, axis=0), test_units, truncation=trunc_rule)
        if np.isfinite(v):
            vals.append(float(v))
    for j in range(n_test):
        v = expected_f2_analytical(times, ref_units, np.delete(test_units, j, axis=0), truncation=trunc_rule)
        if np.isfinite(v):
            vals.append(float(v))
    return np.asarray(vals, dtype=float)


def _bca_ci(
    boot_values: np.ndarray,
    point_est: float,
    jackknife_values: np.ndarray,
    *,
    ci_percentile: float = 90.0,
) -> tuple[float, float]:
    """BCa confidence interval endpoints (Efron 1987)."""
    boot_values = np.asarray(boot_values, dtype=float)
    jackknife_values = np.asarray(jackknife_values, dtype=float)
    n = len(boot_values)
    if n < 2 or len(jackknife_values) < 2:
        return float("nan"), float("nan")

    prop = float(np.mean(boot_values < point_est))
    prop = np.clip(prop, 1.0 / (n + 1), n / (n + 1))
    z0 = float(norm.ppf(prop))

    jmean = float(np.mean(jackknife_values))
    u = jmean - jackknife_values
    denom = 6.0 * (float(np.sum(u**2)) ** 1.5)
    a = float(np.sum(u**3) / denom) if denom > 1e-20 else 0.0

    alpha = (100.0 - ci_percentile) / 100.0

    def _adj_prob(tail: float) -> float:
        z = float(norm.ppf(tail))
        num = z0 + z
        den = 1.0 - a * num
        if abs(den) < 1e-12:
            return tail
        return float(norm.cdf(z0 + num / den))

    lo_p = _adj_prob(alpha / 2.0)
    hi_p = _adj_prob(1.0 - alpha / 2.0)
    lo = float(np.percentile(boot_values, 100.0 * lo_p))
    hi = float(np.percentile(boot_values, 100.0 * hi_p))
    return lo, hi


def bootstrap_ef2_validation(
    times: np.ndarray,
    ref_units: np.ndarray,
    test_units: np.ndarray,
    *,
    n_iterations: int = 10_000,
    seed: int = 42,
    truncation: str = "EMA",
    ci_percentile: float = 90.0,
) -> BootstrapEf2Validation:
    """Unit-level bootstrap of analytical E(f2) with percentile and BCa 90% CIs."""
    ref_units = np.asarray(ref_units, dtype=float)
    test_units = np.asarray(test_units, dtype=float)
    n_ref, n_test = ref_units.shape[0], test_units.shape[0]
    trunc_rule = _resolve_trunc_rule(times, ref_units, test_units, truncation)

    ref_means = np.mean(ref_units, axis=0)
    test_means = np.mean(test_units, axis=0)
    f2_obs = _observed_f2_primary(times, ref_means, test_means)
    point_e = expected_f2_analytical(times, ref_units, test_units, truncation=trunc_rule)

    rng = np.random.default_rng(seed)
    boot_values = np.empty(n_iterations, dtype=float)
    for i in range(n_iterations):
        ref_idx = rng.choice(n_ref, size=n_ref, replace=True)
        test_idx = rng.choice(n_test, size=n_test, replace=True)
        boot_values[i] = expected_f2_analytical(
            times,
            ref_units[ref_idx],
            test_units[test_idx],
            truncation=trunc_rule,
        )

    alpha = (100.0 - ci_percentile) / 2.0
    pct_lo, pct_hi = np.percentile(boot_values, [alpha, 100.0 - alpha])
    jack = _jackknife_ef2_values(times, ref_units, test_units, trunc_rule)
    bca_lo, bca_hi = _bca_ci(boot_values, point_e, jack, ci_percentile=ci_percentile)
    bc = 2.0 * f2_obs - point_e

    return BootstrapEf2Validation(
        f2_observed=float(f2_obs),
        e_f2_point=float(point_e),
        bc_f2=float(bc),
        percentile_ci_lower=float(pct_lo),
        percentile_ci_upper=float(pct_hi),
        bca_ci_lower=float(bca_lo),
        bca_ci_upper=float(bca_hi),
        median=float(np.median(boot_values)),
        n_iterations=n_iterations,
        seed=seed,
        truncation=trunc_rule,
        similar_percentile=float(pct_lo) >= 50.0,
        similar_bca=float(bca_lo) >= 50.0 if np.isfinite(bca_lo) else False,
    )


def bootstrap_ef2_columnwise(
    times: np.ndarray,
    ref_units: np.ndarray,
    test_units: np.ndarray,
    *,
    n_iterations: int = 10_000,
    seed: int = 42,
    truncation: str = "EMA",
    ci_percentile: float = 90.0,
) -> BootstrapEf2Validation:
    """Legacy per-time-point column resampling (wrong granularity) for Table 5 contrast."""
    ref_units = np.asarray(ref_units, dtype=float)
    test_units = np.asarray(test_units, dtype=float)
    n_ref, n_test = ref_units.shape[0], test_units.shape[0]
    n_t = times.shape[0]
    trunc_rule = _resolve_trunc_rule(times, ref_units, test_units, truncation)

    ref_means = np.mean(ref_units, axis=0)
    test_means = np.mean(test_units, axis=0)
    f2_obs = _observed_f2_primary(times, ref_means, test_means)
    point_e = expected_f2_analytical(times, ref_units, test_units, truncation=trunc_rule)

    rng = np.random.default_rng(seed)
    boot_values = np.empty(n_iterations, dtype=float)
    for i in range(n_iterations):
        ref_sample = np.empty((n_ref, n_t))
        test_sample = np.empty((n_test, n_t))
        for t in range(n_t):
            ref_sample[:, t] = rng.choice(ref_units[:, t], size=n_ref, replace=True)
            test_sample[:, t] = rng.choice(test_units[:, t], size=n_test, replace=True)
        boot_values[i] = expected_f2_analytical(
            times, ref_sample, test_sample, truncation=trunc_rule
        )

    alpha = (100.0 - ci_percentile) / 2.0
    pct_lo, pct_hi = np.percentile(boot_values, [alpha, 100.0 - alpha])
    jack = _jackknife_ef2_values(times, ref_units, test_units, trunc_rule)
    bca_lo, bca_hi = _bca_ci(boot_values, point_e, jack, ci_percentile=ci_percentile)
    bc = 2.0 * f2_obs - point_e

    return BootstrapEf2Validation(
        f2_observed=float(f2_obs),
        e_f2_point=float(point_e),
        bc_f2=float(bc),
        percentile_ci_lower=float(pct_lo),
        percentile_ci_upper=float(pct_hi),
        bca_ci_lower=float(bca_lo),
        bca_ci_upper=float(bca_hi),
        median=float(np.median(boot_values)),
        n_iterations=n_iterations,
        seed=seed,
        truncation=trunc_rule,
        similar_percentile=float(pct_lo) >= 50.0,
        similar_bca=float(bca_lo) >= 50.0 if np.isfinite(bca_lo) else False,
    )
