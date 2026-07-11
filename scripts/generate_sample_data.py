"""Generate sample dissolution CSV datasets for the reference benchmark."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd


def weibull(t: np.ndarray, alpha: float, beta: float) -> np.ndarray:
    return 100.0 * (1.0 - np.exp(-np.power(t, beta) / alpha))


def write_csv(path: Path, times: np.ndarray, units: np.ndarray, n_reps: int) -> None:
    cols = {"Time": times}
    for i in range(n_reps):
        cols[f"Rep {i + 1}"] = units[i]
    pd.DataFrame(cols).to_csv(path, index=False)


def main() -> None:
    data_dir = Path(__file__).resolve().parent.parent / "data"
    data_dir.mkdir(parents=True, exist_ok=True)

    # Dense early sampling supports the <=60% K-P window (>=3 points per unit).
    times = np.array([5.0, 10.0, 15.0, 20.0, 25.0, 30.0, 45.0, 60.0, 90.0, 120.0])
    rng = np.random.default_rng(2024)

    # Low-variability reference: target Weibull alpha~58.97, beta~1.4865, R2~0.996
    alpha_ref, beta_ref = 58.97, 1.4865
    base_ref = weibull(times, alpha_ref, beta_ref)
    n_ref = 6
    noise_ref = rng.normal(0, 2.2, size=(n_ref, len(times)))
    units_ref = np.clip(base_ref + noise_ref, 0, 100)
    write_csv(data_dir / "weibull_ref.csv", times, units_ref, n_ref)

    # Low-variability test: near-identical profile for f2~91
    alpha_test, beta_test = 59.35, 1.4840
    base_test = weibull(times, alpha_test, beta_test)
    n_test = 6
    noise_test = rng.normal(0, 2.2, size=(n_test, len(times)))
    units_test = np.clip(base_test + noise_test, 0, 100)
    write_csv(data_dir / "weibull_test.csv", times, units_test, n_test)

    # High-variability pair for bootstrap gate testing
    alpha_hi, beta_hi = 55.0, 1.5
    base_hi = weibull(times, alpha_hi, beta_hi)
    n_hi = 8
    noise_hi = rng.normal(0, 5.5, size=(n_hi, len(times)))
    units_hi_ref = np.clip(base_hi + noise_hi, 0, 100)
    units_hi_test = np.clip(base_hi + noise_hi * 1.15 + 2.0, 0, 100)
    write_csv(data_dir / "high_cv_ref.csv", times, units_hi_ref, n_hi)
    write_csv(data_dir / "high_cv_test.csv", times, units_hi_test, n_hi)

    print(f"Sample datasets written to {data_dir}")


if __name__ == "__main__":
    main()
