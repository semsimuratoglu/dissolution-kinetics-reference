"""Similarity, difference, and exposure metrics for dissolution profiles."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class TruncationResult:
    """Result of a profile truncation rule."""

    times: np.ndarray
    ref_means: np.ndarray
    test_means: np.ndarray
    kept_indices: np.ndarray
    rule: str


def fda_truncation(
    times: np.ndarray, ref_means: np.ndarray, test_means: np.ndarray
) -> TruncationResult:
    """FDA pathway: keep points up to and including the first where BOTH means exceed 85%."""
    both_above = (ref_means > 85.0) & (test_means > 85.0)
    if np.any(both_above):
        last_idx = int(np.argmax(both_above))
        kept = np.arange(last_idx + 1)
    else:
        kept = np.arange(len(times))
    return TruncationResult(
        times=times[kept],
        ref_means=ref_means[kept],
        test_means=test_means[kept],
        kept_indices=kept,
        rule="FDA",
    )


def ema_truncation(
    times: np.ndarray, ref_means: np.ndarray, test_means: np.ndarray
) -> TruncationResult:
    """EMA pathway: keep points up to and including the first where EITHER mean exceeds 85%."""
    any_above = (ref_means > 85.0) | (test_means > 85.0)
    if np.any(any_above):
        last_idx = int(np.argmax(any_above))
        kept = np.arange(last_idx + 1)
    else:
        kept = np.arange(len(times))
    return TruncationResult(
        times=times[kept],
        ref_means=ref_means[kept],
        test_means=test_means[kept],
        kept_indices=kept,
        rule="EMA",
    )


def f2(ref_means: np.ndarray, test_means: np.ndarray) -> float:
    """f2 similarity factor.

    f2 = 50 * log10( 100 / sqrt(1 + (1/L) * sum((Rbar_t - Tbar_t)**2)) )
    """
    ref = np.asarray(ref_means, dtype=float)
    test = np.asarray(test_means, dtype=float)
    if ref.shape != test.shape:
        raise ValueError("ref_means and test_means must have the same length")
    length = len(ref)
    if length == 0:
        return float("nan")
    mean_sq_diff = float(np.mean((ref - test) ** 2))
    inner = np.sqrt(1.0 + mean_sq_diff)
    return float(50.0 * np.log10(100.0 / inner))


def expected_f2(ref_units: np.ndarray, test_units: np.ndarray) -> float:
    """Expected f2 E(f2) per Xu 2021 Eq.10 / Noce 2020 Eq.3.

    E(f2) = 100 - 25*log10(1 + (1/P)*[ sum((Rbar-Tbar)^2) + (1/n)*sum(s2_R + s2_T) ])
    """
    ref_units = np.asarray(ref_units, dtype=float)
    test_units = np.asarray(test_units, dtype=float)
    ref_means = np.mean(ref_units, axis=0)
    test_means = np.mean(test_units, axis=0)
    p = ref_means.shape[0]
    if p == 0:
        return float("nan")
    n_ref, n_test = ref_units.shape[0], test_units.shape[0]
    sum_sq = float(np.sum((ref_means - test_means) ** 2))
    var_term = 0.0
    for j in range(p):
        if n_ref > 1:
            var_term += float(np.var(ref_units[:, j], ddof=1)) / n_ref
        if n_test > 1:
            var_term += float(np.var(test_units[:, j], ddof=1)) / n_test
    inner = sum_sq / p + var_term
    return float(100.0 - 25.0 * np.log10(1.0 + inner))


def f1(ref_means: np.ndarray, test_means: np.ndarray) -> float:
    """f1 difference factor on the same retained slice as f2.

    f1 = 100 * sum(|Rbar_t - Tbar_t|) / sum(Rbar_t)
    """
    ref = np.asarray(ref_means, dtype=float)
    test = np.asarray(test_means, dtype=float)
    ref_sum = float(np.sum(ref))
    if ref_sum == 0:
        return float("nan")
    return float(100.0 * np.sum(np.abs(ref - test)) / ref_sum)


def de_percent(times: np.ndarray, means: np.ndarray) -> float:
    """DE% on the full measured profile.

    Prepends (t0=0, Q0=0), integrates with ``numpy.trapezoid``, and divides by
    ``100 * t_last`` where ``t_last`` is the last measured time point.
    """
    times = np.asarray(times, dtype=float)
    means = np.asarray(means, dtype=float)
    if len(times) == 0:
        return float("nan")
    t_ext = np.concatenate(([0.0], times))
    q_ext = np.concatenate(([0.0], means))
    auc = float(np.trapezoid(q_ext, t_ext))
    t_last = float(times[-1])
    if t_last == 0:
        return float("nan")
    return float(auc / (100.0 * t_last) * 100.0)


def mdt(times: np.ndarray, means: np.ndarray) -> float:
    """MDT on the full measured profile (DDSolver midpoint-increment convention).

    With Q0=0 at t0=0 prepended:
      dQ_i = Q(t_i) - Q(t_{i-1}),  tmid_i = (t_i + t_{i-1}) / 2
      MDT = sum(tmid_i * dQ_i) / sum(dQ_i)
    """
    times = np.asarray(times, dtype=float)
    means = np.asarray(means, dtype=float)
    if len(times) == 0:
        return float("nan")
    t_ext = np.concatenate(([0.0], times))
    q_ext = np.concatenate(([0.0], means))
    dQ = np.diff(q_ext)
    t_mid = 0.5 * (t_ext[:-1] + t_ext[1:])
    dQ_sum = float(np.sum(dQ))
    if dQ_sum == 0:
        return float("nan")
    return float(np.sum(t_mid * dQ) / dQ_sum)


def expected_f2_analytical(
    times: np.ndarray,
    ref_units: np.ndarray,
    test_units: np.ndarray,
    *,
    truncation: str | None = None,
) -> float:
    """Expected f2 estimator E(f2) per Xu 2021 Eq. 10 / Noce 2020 Eq. 3.

    Uses EMA-primary truncation when FDA and EMA slices diverge.
  """
    times = np.asarray(times, dtype=float)
    ref_units = np.asarray(ref_units, dtype=float)
    test_units = np.asarray(test_units, dtype=float)
    n_ref, n_test = ref_units.shape[0], test_units.shape[0]
    if ref_units.shape[1] < 3:
        return float("nan")

    ref_means = np.mean(ref_units, axis=0)
    test_means = np.mean(test_units, axis=0)
    fda = fda_truncation(times, ref_means, test_means)
    ema = ema_truncation(times, ref_means, test_means)
    slices_diverge = not np.array_equal(fda.kept_indices, ema.kept_indices)
    if truncation is not None:
        trunc = ema if truncation.upper() == "EMA" else fda
    elif slices_diverge:
        trunc = ema
    else:
        trunc = fda

    kept = trunc.kept_indices
    p = len(kept)
    if p < 3:
        return float("nan")

    ref_t = trunc.ref_means
    test_t = trunc.test_means
    diff_sq_sum = float(np.sum((ref_t - test_t) ** 2))
    var_sum = 0.0
    for j in kept:
        if n_ref > 1:
            s_ref = float(np.std(ref_units[:, j], ddof=1))
            var_sum += (s_ref * s_ref) / n_ref
        if n_test > 1:
            s_test = float(np.std(test_units[:, j], ddof=1))
            var_sum += (s_test * s_test) / n_test

    inner = (1.0 / p) * (diff_sq_sum + var_sum)
    if not np.isfinite(inner) or inner < 0:
        return float("nan")
    return float(100.0 - 25.0 * np.log10(1.0 + inner))


def compute_profile_metrics(
    times: np.ndarray,
    ref_units: np.ndarray,
    test_units: np.ndarray,
    truncation: str = "FDA",
) -> dict:
    """Compute f1, f2, DE%, and MDT for a reference/test pair.

    f1 and f2 use the FDA or EMA 85% truncated slice. DE% and MDT always use the
    full measured mean profile (pathway-independent).
    """
    ref_means = np.mean(ref_units, axis=0)
    test_means = np.mean(test_units, axis=0)

    if truncation.upper() == "EMA":
        trunc = ema_truncation(times, ref_means, test_means)
    else:
        trunc = fda_truncation(times, ref_means, test_means)

    return {
        "truncation": trunc.rule,
        "kept_indices": trunc.kept_indices.tolist(),
        "f2": f2(trunc.ref_means, trunc.test_means),
        "f1": f1(trunc.ref_means, trunc.test_means),
        "DE_ref": de_percent(times, ref_means),
        "DE_test": de_percent(times, test_means),
        "MDT_ref": mdt(times, ref_means),
        "MDT_test": mdt(times, test_means),
        "truncated_times": trunc.times.tolist(),
        "ref_means_truncated": trunc.ref_means.tolist(),
        "test_means_truncated": trunc.test_means.tolist(),
    }
