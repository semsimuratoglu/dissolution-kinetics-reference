"""Unit-level bootstrap resampling for expected f2 (E-f2) similarity uncertainty."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from refeng.metrics import expected_f2_analytical, f2, fda_truncation, ema_truncation


@dataclass
class BootstrapF2Result:
    """Bootstrap summary for E(f2)."""

    f2_observed: float
    """Point analytical E(f2) on observed unit data (Xu 2021 / Noce 2020)."""
    e_f2: float
    """Bootstrap mean of resampled E(f2) values (diagnostic)."""
    mean: float
    median: float
    ci_lower: float
    ci_upper: float
    bc_f2: float
    similar: bool
    n_iterations: int
    truncation: str
    seed: int


def _primary_truncation_rule(times: np.ndarray, ref_means: np.ndarray, test_means: np.ndarray) -> str:
    fda = fda_truncation(times, ref_means, test_means)
    ema = ema_truncation(times, ref_means, test_means)
    if np.array_equal(fda.kept_indices, ema.kept_indices):
        return "FDA"
    return "EMA"


def _observed_f2_primary(times: np.ndarray, ref_means: np.ndarray, test_means: np.ndarray) -> float:
    rule = _primary_truncation_rule(times, ref_means, test_means)
    if rule == "EMA":
        trunc = ema_truncation(times, ref_means, test_means)
    else:
        trunc = fda_truncation(times, ref_means, test_means)
    return f2(trunc.ref_means, trunc.test_means)


def bootstrap_f2(
    times: np.ndarray,
    ref_units: np.ndarray,
    test_units: np.ndarray,
    *,
    n_iterations: int = 10_000,
    seed: int = 42,
    truncation: str = "FDA",
    ci_percentile: float = 90.0,
) -> BootstrapF2Result:
    """Unit-level bootstrap for E(f2) with bias-corrected point estimate.

    Each iteration resamples whole experimental units (rows) with replacement,
    independently for reference and test. The 90% percentile CI is taken from the
    bootstrap distribution of analytical E(f2) (not observed f2).

    Similarity decision: lower CI bound of E(f2) >= 50.
    bc-f2 = 2 * f2_observed - E(f2) on observed data.
    """
    ref_units = np.asarray(ref_units, dtype=float)
    test_units = np.asarray(test_units, dtype=float)
    n_ref, n_test = ref_units.shape[0], test_units.shape[0]

    ref_means_obs = np.mean(ref_units, axis=0)
    test_means_obs = np.mean(test_units, axis=0)
    rule = _primary_truncation_rule(times, ref_means_obs, test_means_obs)
    if truncation.upper() == "EMA" or rule == "EMA":
        trunc_rule = "EMA"
    else:
        trunc_rule = "FDA"

    f2_obs = _observed_f2_primary(times, ref_means_obs, test_means_obs)
    point_e = expected_f2_analytical(times, ref_units, test_units, truncation=trunc_rule)

    rng = np.random.default_rng(seed)
    boot_values = np.empty(n_iterations, dtype=float)

    ref_indices_orig = np.arange(n_ref)
    test_indices_orig = np.arange(n_test)

    for i in range(n_iterations):
        ref_idx = rng.choice(ref_indices_orig, size=n_ref, replace=True)
        test_idx = rng.choice(test_indices_orig, size=n_test, replace=True)
        boot_values[i] = expected_f2_analytical(
            times,
            ref_units[ref_idx],
            test_units[test_idx],
            truncation=trunc_rule,
        )

    alpha = (100.0 - ci_percentile) / 2.0
    ci_lo, ci_hi = np.percentile(boot_values, [alpha, 100.0 - alpha])
    boot_mean = float(np.mean(boot_values))
    boot_median = float(np.median(boot_values))
    bc = 2.0 * f2_obs - point_e

    return BootstrapF2Result(
        f2_observed=float(f2_obs),
        e_f2=float(point_e),
        mean=boot_mean,
        median=boot_median,
        ci_lower=float(ci_lo),
        ci_upper=float(ci_hi),
        bc_f2=float(bc),
        similar=float(ci_lo) >= 50.0,
        n_iterations=n_iterations,
        truncation=trunc_rule,
        seed=seed,
    )


def resample_unit_indices(
    n_units: int,
    rng: np.random.Generator,
) -> np.ndarray:
    """Draw a unit-level bootstrap sample (test hook)."""
    return rng.choice(np.arange(n_units), size=n_units, replace=True)
