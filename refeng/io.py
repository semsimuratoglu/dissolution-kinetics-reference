"""CSV loading utilities for dissolution unit matrices."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd


def load_dissolution_csv(path: str | Path) -> tuple[np.ndarray, np.ndarray]:
    """Load a dissolution CSV into time vector and unit matrix.

    Expected columns: ``Time`` followed by replicate columns (``Rep 1``, ``Rep 2``, …).

    Parameters
    ----------
    path:
        Path to the CSV file.

    Returns
    -------
    times : ndarray, shape (T,)
        Sorted, unique time points in ascending order.
    units : ndarray, shape (N, T)
        Percent dissolved for each experimental unit (row) at each time point (column).
        Missing values are not supported; all replicate columns must be numeric.
    """
    path = Path(path)
    df = pd.read_csv(path)

    time_col = None
    for col in df.columns:
        if str(col).strip().lower().startswith("time"):
            time_col = col
            break
    if time_col is None:
        raise ValueError(f"{path}: expected a 'Time' column")
    if time_col != "Time":
        df = df.rename(columns={time_col: "Time"})

    rep_cols = [c for c in df.columns if c != "Time"]
    if not rep_cols:
        raise ValueError(f"{path}: expected at least one replicate column")

    times = df["Time"].to_numpy(dtype=float)
    values = df[rep_cols].to_numpy(dtype=float)

    order = np.argsort(times)
    times = times[order]
    values = values[order]

    # Each replicate column is one experimental unit; transpose to (N_units, T).
    units = values.T
    return times, units


def unit_means(units: np.ndarray) -> np.ndarray:
    """Return the mean profile across units at each time point."""
    return np.mean(units, axis=0)


def coefficient_of_variation(units: np.ndarray) -> np.ndarray:
    """Per-time-point %CV across units."""
    n_units = units.shape[0]
    if n_units < 2:
        return np.full(units.shape[1], np.nan)
    means = np.mean(units, axis=0)
    stds = np.std(units, axis=0, ddof=1)
    with np.errstate(divide="ignore", invalid="ignore"):
        cv = np.where(means != 0, 100.0 * stds / np.abs(means), np.nan)
    return cv


def cv_gate_exceeded(units: np.ndarray) -> tuple[bool, str]:
    """Check whether bootstrap is required by the CV gate.

    First time point: %CV > 20%.
    All later time points: %CV > 10%.
    """
    cv = coefficient_of_variation(units)
    if len(cv) == 0:
        return False, ""

    flags: list[str] = []
    if not np.isnan(cv[0]) and cv[0] > 20.0:
        flags.append(f"t0 CV={cv[0]:.2f}% > 20%")
    for i in range(1, len(cv)):
        if not np.isnan(cv[i]) and cv[i] > 10.0:
            flags.append(f"t{i} CV={cv[i]:.2f}% > 10%")
    return bool(flags), "; ".join(flags)
