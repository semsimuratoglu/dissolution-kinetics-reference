"""Kinetic model equations and per-unit nonlinear least-squares fitting."""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy.optimize import curve_fit

from refeng.gof import goodness_of_fit


@dataclass
class UnitFitResult:
    """Fit result for a single experimental unit."""

    unit_index: int
    params: dict[str, float]
    converged: bool
    message: str = ""
    gof: dict | None = None


@dataclass
class ModelFitSummary:
    """Aggregated fit across all units for one kinetic model."""

    model_name: str
    unit_results: list[UnitFitResult] = field(default_factory=list)
    param_mean: dict[str, float] = field(default_factory=dict)
    param_sd: dict[str, float] = field(default_factory=dict)
    n_units_total: int = 0
    n_units_converged: int = 0
    gof_mean: dict[str, float | None] = field(default_factory=dict)
    notes: str = ""


# --- Model equations (F in percent dissolved, t in minutes) ---


def zero_order(t: np.ndarray, k0: float) -> np.ndarray:
    """Zero-order: F = k0 * t."""
    return k0 * t


def first_order(t: np.ndarray, k1: float) -> np.ndarray:
    """First-order: F = 100 * (1 - exp(-k1 * t))."""
    return 100.0 * (1.0 - np.exp(-k1 * t))


def higuchi(t: np.ndarray, kH: float) -> np.ndarray:
    """Higuchi: F = kH * sqrt(t)."""
    return kH * np.sqrt(np.maximum(t, 0.0))


def hixson_crowell(t: np.ndarray, kHC: float) -> np.ndarray:
    """Hixson-Crowell: F = 100 * (1 - (1 - kHC * t)**3), unconstrained NLS on all points."""
    return 100.0 * (1.0 - (1.0 - kHC * t) ** 3)


def korsmeyer_peppas(t: np.ndarray, kKP: float, n: float) -> np.ndarray:
    """Korsmeyer-Peppas: F = kKP * t**n."""
    return kKP * np.power(np.maximum(t, 0.0), n)


def weibull(t: np.ndarray, alpha: float, beta: float) -> np.ndarray:
    """Weibull (scale form): F = 100 * (1 - exp(-(t**beta) / alpha))."""
    tb = np.power(np.maximum(t, 0.0), beta)
    return 100.0 * (1.0 - np.exp(-tb / alpha))


MODEL_SPECS: dict[str, dict] = {
    "zero_order": {
        "func": zero_order,
        "param_names": ["k0"],
        "p0": [1.0],
        "bounds": ([0.0], [np.inf]),
        "positive": True,
    },
    "first_order": {
        "func": first_order,
        "param_names": ["k1"],
        "p0": [0.05],
        "bounds": ([1e-8], [np.inf]),
        "positive": True,
    },
    "higuchi": {
        "func": higuchi,
        "param_names": ["kH"],
        "p0": [5.0],
        "bounds": ([0.0], [np.inf]),
        "positive": True,
    },
    "hixson_crowell": {
        "func": hixson_crowell,
        "param_names": ["kHC"],
        "p0": [0.01],
        "bounds": ([1e-8], [1.0]),
        "positive": True,
    },
    "korsmeyer_peppas_full": {
        "func": korsmeyer_peppas,
        "param_names": ["kKP", "n"],
        "p0": [10.0, 0.5],
        "bounds": ([1e-8, 0.0], [np.inf, 10.0]),
        "positive": True,
    },
    "korsmeyer_peppas_60": {
        "func": korsmeyer_peppas,
        "param_names": ["kKP", "n"],
        "p0": [10.0, 0.5],
        "bounds": ([1e-8, 0.0], [np.inf, 10.0]),
        "positive": True,
        "max_release": 60.0,
    },
    "weibull": {
        "func": weibull,
        "param_names": ["alpha", "beta"],
        "p0": [50.0, 1.5],
        "bounds": ([1e-8, 1e-8], [np.inf, 10.0]),
        "positive": True,
    },
}


def _fit_single_unit(
    times: np.ndarray,
    observed: np.ndarray,
    model_key: str,
) -> UnitFitResult:
    """Fit one kinetic model to one unit's dissolution profile."""
    spec = MODEL_SPECS[model_key]
    func = spec["func"]
    param_names: list[str] = spec["param_names"]
    n_params = len(param_names)

    t_fit = times
    y_fit = observed

    if "max_release" in spec:
        mask = observed <= spec["max_release"]
        if np.sum(mask) < 3:
            return UnitFitResult(
                unit_index=-1,
                params={},
                converged=False,
                message=f"< 3 points with release <= {spec['max_release']}%",
            )
        t_fit = times[mask]
        y_fit = observed[mask]

    p0 = spec["p0"]
    bounds = spec["bounds"]
    positive = spec.get("positive", False)

    try:
        if positive:
            popt, _pcov = curve_fit(
                func,
                t_fit,
                y_fit,
                p0=p0,
                bounds=bounds,
                method="trf",
                maxfev=10000,
            )
        else:
            popt, _pcov = curve_fit(
                func,
                t_fit,
                y_fit,
                p0=p0,
                method="lm",
                maxfev=10000,
            )
        params = {name: float(val) for name, val in zip(param_names, popt)}
        predicted = func(t_fit, *popt)
        gof = goodness_of_fit(y_fit, predicted, n_params)
        return UnitFitResult(
            unit_index=-1,
            params=params,
            converged=True,
            gof=gof,
        )
    except Exception as exc:
        return UnitFitResult(
            unit_index=-1,
            params={},
            converged=False,
            message=str(exc),
        )


def fit_model_per_unit(
    times: np.ndarray,
    units: np.ndarray,
    model_key: str,
) -> ModelFitSummary:
    """Fit a kinetic model independently to each unit, then aggregate."""
    spec = MODEL_SPECS[model_key]
    param_names: list[str] = spec["param_names"]
    n_units = units.shape[0]

    unit_results: list[UnitFitResult] = []
    for i in range(n_units):
        result = _fit_single_unit(times, units[i], model_key)
        result.unit_index = i
        unit_results.append(result)

    converged = [r for r in unit_results if r.converged]
    n_conv = len(converged)

    summary = ModelFitSummary(
        model_name=model_key,
        unit_results=unit_results,
        n_units_total=n_units,
        n_units_converged=n_conv,
    )

    if model_key == "korsmeyer_peppas_60":
        summary.notes = f"converged {n_conv}/{n_units} units (<=60% window)"

    if n_conv == 0:
        return summary

    for name in param_names:
        vals = [r.params[name] for r in converged]
        summary.param_mean[name] = float(np.mean(vals))
        summary.param_sd[name] = float(np.std(vals, ddof=1)) if len(vals) > 1 else 0.0

    gof_keys = ["R2", "adj_R2", "RMSE", "AIC", "MSC"]
    for key in gof_keys:
        vals = [r.gof[key] for r in converged if r.gof and r.gof[key] is not None]
        if vals:
            summary.gof_mean[key] = float(np.mean(vals))
        else:
            summary.gof_mean[key] = None

    return summary


def fit_all_models(times: np.ndarray, units: np.ndarray) -> dict[str, ModelFitSummary]:
    """Fit every registered kinetic model to the unit matrix."""
    return {
        key: fit_model_per_unit(times, units, key)
        for key in MODEL_SPECS
    }
