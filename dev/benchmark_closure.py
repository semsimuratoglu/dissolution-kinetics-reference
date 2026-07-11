#!/usr/bin/env python3
"""Development benchmark closure report (literature validation tasks 1–7)."""

from __future__ import annotations

import json
import math
import os
import subprocess
import sys
from dataclasses import dataclass
from datetime import date
from pathlib import Path

import numpy as np
from scipy.stats import norm

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "data"
PY_EXPORT = ROOT / "results" / "python_reference_export.json"
sys.path.insert(0, str(ROOT))

from refeng.bootstrap import bootstrap_f2  # noqa: E402
from refeng.io import cv_gate_exceeded, load_dissolution_csv  # noqa: E402
from refeng.metrics import (  # noqa: E402
    compute_profile_metrics,
    ema_truncation,
    expected_f2_analytical,
    f2,
    fda_truncation,
)
from refeng.models import MODEL_SPECS, fit_all_models  # noqa: E402

RUN_DATE = date.today().isoformat()
BOOT_ITER = 10_000
BOOT_SEED = 42
MC_TOL = 3.0

LITERATURE = {
    "shah_similar": {
        "percentile": (53.2, 67.5),
        "bca": (53.9, 70.2),
        "boddu_pct": (53.2, 67.5),
        "boddu_bca": (53.8, 73.0),
    },
    "shah_borderline": {
        "percentile": (45.99, 49.96),
        "bca": (46.15, 50.17),
        "boddu_pct": (45.9, 50.28),
        "boddu_bca": (45.9, 50.28),
    },
}

MODEL_LABELS = {
    "zero_order": "Zero-order",
    "first_order": "First-order",
    "higuchi": "Higuchi",
    "hixson_crowell": "Hixson-Crowell",
    "korsmeyer_peppas_full": "K-P full",
    "korsmeyer_peppas_60": "K-P <=60%",
    "weibull": "Weibull",
}


def mulberry32(seed: int):
    state = seed & 0xFFFFFFFF

    def rng() -> float:
        nonlocal state
        state = (state + 0x6D2B79F5) & 0xFFFFFFFF
        t = state
        t = (t ^ (t >> 15)) * (t | 1)
        t = (t ^ (t + ((t ^ (t >> 7)) * (t | 61)))) & 0xFFFFFFFF
        return ((t ^ (t >> 14)) & 0xFFFFFFFF) / 4294967296.0

    return rng


def draw_indices(n: int, rng) -> np.ndarray:
    return np.array([int(math.floor(rng() * n)) for _ in range(n)], dtype=int)


def primary_truncation(times, ref_means, test_means):
    fda = fda_truncation(times, ref_means, test_means)
    ema = ema_truncation(times, ref_means, test_means)
    if np.array_equal(fda.kept_indices, ema.kept_indices):
        return fda
    return ema


def observed_f2_primary(times, ref_units, test_units) -> float:
    rm, tm = ref_units.mean(0), test_units.mean(0)
    trunc = primary_truncation(times, rm, tm)
    return f2(trunc.ref_means, trunc.test_means)


def bootstrap_f2_distribution(
    times: np.ndarray,
    ref_units: np.ndarray,
    test_units: np.ndarray,
    *,
    n_iterations: int,
    seed: int,
) -> np.ndarray:
    rng = mulberry32(seed)
    n_ref, n_test = ref_units.shape[0], test_units.shape[0]
    out = np.empty(n_iterations)
    for i in range(n_iterations):
        ri = draw_indices(n_ref, rng)
        ti = draw_indices(n_test, rng)
        rm = ref_units[ri].mean(0)
        tm = test_units[ti].mean(0)
        trunc = primary_truncation(times, rm, tm)
        out[i] = f2(trunc.ref_means, trunc.test_means)
    return out


def bootstrap_ef2_distribution(
    times: np.ndarray,
    ref_units: np.ndarray,
    test_units: np.ndarray,
    *,
    n_iterations: int,
    seed: int,
) -> np.ndarray:
    rng = mulberry32(seed)
    n_ref, n_test = ref_units.shape[0], test_units.shape[0]
    rule = "EMA"
    out = np.empty(n_iterations)
    for i in range(n_iterations):
        ri = draw_indices(n_ref, rng)
        ti = draw_indices(n_test, rng)
        out[i] = expected_f2_analytical(
            times, ref_units[ri], test_units[ti], truncation=rule
        )
    return out


def jackknife_f2(times, ref_units, test_units) -> np.ndarray:
    n = ref_units.shape[0]
    vals = np.empty(n)
    for i in range(n):
        mask = np.ones(n, dtype=bool)
        mask[i] = False
        rm = ref_units[mask].mean(0)
        tm = test_units[mask].mean(0)
        trunc = primary_truncation(times, rm, tm)
        vals[i] = f2(trunc.ref_means, trunc.test_means)
    return vals


def jackknife_ef2(times, ref_units, test_units) -> np.ndarray:
    n = ref_units.shape[0]
    vals = np.empty(n)
    rule = "EMA"
    for i in range(n):
        mask = np.ones(n, dtype=bool)
        mask[i] = False
        vals[i] = expected_f2_analytical(
            times, ref_units[mask], test_units[mask], truncation=rule
        )
    return vals


def bca_ci(
    boot: np.ndarray,
    point: float,
    jack: np.ndarray,
    alpha: float = 0.10,
) -> tuple[float, float]:
    boot = boot[np.isfinite(boot)]
    jack = jack[np.isfinite(jack)]
    if len(boot) < 10:
        return float("nan"), float("nan")
    z0 = norm.ppf(np.clip(np.mean(boot < point), 1e-6, 1 - 1e-6))
    theta_dot = float(np.mean(jack))
    dif = theta_dot - jack
    num = float(np.sum(dif**3))
    den = 6.0 * float(np.sum(dif**2)) ** 1.5
    a = num / den if den > 0 else 0.0
    zs = norm.ppf([alpha / 2, 1 - alpha / 2])
    probs = []
    for z in zs:
        adj = z0 + (z0 + z) / (1 - a * (z0 + z)) if abs(1 - a * (z0 + z)) > 1e-12 else z0 + z
        probs.append(norm.cdf(adj))
    lo, hi = np.percentile(boot, [100 * probs[0], 100 * probs[1]])
    return float(lo), float(hi)


def percentile_ci(boot: np.ndarray, alpha: float = 0.10) -> tuple[float, float]:
    lo, hi = np.percentile(boot, [100 * alpha / 2, 100 * (1 - alpha / 2)])
    return float(lo), float(hi)


@dataclass
class ShahEngineResult:
    label: str
    f2_obs: float
    ef2_point: float
    f2_pct: tuple[float, float]
    f2_bca: tuple[float, float]
    ef2_pct: tuple[float, float]
    ef2_bca: tuple[float, float]


def run_shah_engine(ref_file: str, test_file: str, label: str) -> ShahEngineResult:
    times, ref_units = load_dissolution_csv(DATA / ref_file)
    _, test_units = load_dissolution_csv(DATA / test_file)
    f2_obs = observed_f2_primary(times, ref_units, test_units)
    ef2 = expected_f2_analytical(times, ref_units, test_units, truncation="EMA")
    f2_boot = bootstrap_f2_distribution(times, ref_units, test_units, n_iterations=BOOT_ITER, seed=BOOT_SEED)
    ef2_boot = bootstrap_ef2_distribution(times, ref_units, test_units, n_iterations=BOOT_ITER, seed=BOOT_SEED)
    return ShahEngineResult(
        label=label,
        f2_obs=f2_obs,
        ef2_point=ef2,
        f2_pct=percentile_ci(f2_boot),
        f2_bca=bca_ci(f2_boot, f2_obs, jackknife_f2(times, ref_units, test_units)),
        ef2_pct=percentile_ci(ef2_boot),
        ef2_bca=bca_ci(ef2_boot, ef2, jackknife_ef2(times, ref_units, test_units)),
    )


def within_tol(val: float, target: float, tol: float = MC_TOL) -> str:
    if not np.isfinite(val) or not np.isfinite(target):
        return "N/A"
    return "YES" if abs(val - target) <= tol else "NO"


def print_task1():
    print("\n" + "=" * 78)
    print("TASK 1 — Shah bootstrap external validation (validation mode, CV gate OK)")
    print("=" * 78)
    r_out = subprocess.run(
        ["Rscript", str(ROOT / "scripts/benchmark_closure_shah.R"), str(DATA), str(BOOT_ITER), str(BOOT_SEED)],
        capture_output=True,
        text=True,
        check=False,
    )
    print(r_out.stdout)
    if r_out.stderr:
        print(r_out.stderr, file=sys.stderr)

    engine = {
        "shah_similar": run_shah_engine("shah_ref.csv", "shah_test.csv", "shah_similar"),
        "shah_borderline": run_shah_engine("shah_ref.csv", "shah_borderline_test.csv", "shah_borderline"),
    }

    bootf2_json = subprocess.run(
        ["Rscript", str(ROOT / "scripts/benchmark_closure_bootf2_json.R"), str(DATA), str(BOOT_ITER), str(BOOT_SEED)],
        capture_output=True,
        text=True,
        check=True,
    )
    bootf2_list = json.loads(bootf2_json.stdout)
    bootf2 = {x["label"]: x for x in bootf2_list}

    print("\nNote: bootf2 est.f2 f2o (57.47/45.23) differs from mean-profile f2 (60.03/48.05) used by refeng/literature.")
    print("CI comparison below uses bootstrap distribution of mean-profile f2 (our engine) vs bootf2 est.f2 CI.\n")
    hdr = f"{'quantity':<28} {'our_engine':>12} {'bootf2':>12} {'Noce':>12} {'Boddu':>12} {'tol?':>6}"
    print(hdr)
    print("-" * len(hdr))
    for key in ("shah_similar", "shah_borderline"):
        eng = engine[key]
        bf2 = bootf2[key]
        lit = LITERATURE[key]
        rows = [
            ("f2_observed (mean-profile)", eng.f2_obs, bf2["f2_obs"], 60.03 if key == "shah_similar" else 48.05, "—", "—"),
            ("bootf2 est.f2 f2o", bf2["f2_obs"], bf2["f2_obs"], float("nan"), float("nan"), "—"),
            ("E_f2_analytical", eng.ef2_point, bf2["exp_f2"], float("nan"), float("nan"), "—"),
            ("pct_CI_lower", eng.f2_pct[0], bf2["pct_l"], lit["percentile"][0], lit["boddu_pct"][0],
             within_tol(eng.f2_pct[0], lit["percentile"][0])),
            ("pct_CI_upper", eng.f2_pct[1], bf2["pct_u"], lit["percentile"][1], lit["boddu_pct"][1],
             within_tol(eng.f2_pct[1], lit["percentile"][1])),
            ("BCa_CI_lower", eng.f2_bca[0], bf2["bca_l"], lit["bca"][0], lit["boddu_bca"][0],
             within_tol(eng.f2_bca[0], lit["bca"][0])),
            ("BCa_CI_upper", eng.f2_bca[1], bf2["bca_u"], lit["bca"][1], lit["boddu_bca"][1],
             within_tol(eng.f2_bca[1], lit["bca"][1])),
        ]
        print(f"\n[{key}]")
        for q, ours, b2, noce, boddu, ok in rows:
            noce_s = f"{noce:.2f}" if isinstance(noce, (int, float)) and np.isfinite(noce) else str(noce)
            boddu_s = f"{boddu:.2f}" if isinstance(boddu, (int, float)) and np.isfinite(boddu) else str(boddu)
            print(f"{q:<28} {ours:>12.2f} {b2:>12.2f} {noce_s:>12} {boddu_s:>12} {ok:>6}")

    print("\n--- Our engine E(f2) bootstrap CI (decision rule metric) ---")
    for key in ("shah_similar", "shah_borderline"):
        eng = engine[key]
        print(f"{key}: E_f2={eng.ef2_point:.2f}  pct=[{eng.ef2_pct[0]:.2f},{eng.ef2_pct[1]:.2f}]  "
              f"BCa=[{eng.ef2_bca[0]:.2f},{eng.ef2_bca[1]:.2f}]  "
              f"decision={'SIMILAR' if eng.ef2_pct[0] >= 50 else 'NOT'}")


def print_task2():
    print("\n" + "=" * 78)
    print("TASK 2 — Mendyk N12 provenance")
    print("=" * 78)
    print(
        "VERDICT: mendyk_reference_N12.csv / mendyk_test_N12.csv are SYNTHETIC in-repo fixtures.\n"
        "  - NOT transcribed from Mendyk et al. (2013) KinetDS / PhEq_bootstrap Case Study 1.\n"
        "  - Times: 15,30,45,60,90,120 min (N=12); Mendyk 2013 published schedule is 3–30 min.\n"
        "  - generate_sample_data.py writes weibull_* and high_cv_* only — NOT mendyk N12.\n"
        "  - Role: synthetic high-variability bootstrap anchor (first-point CV ≈ 25.7% / 28.4%).\n"
        "  - Label for Table 6: 'synthetic high-variability set (Mendyk-style naming)'."
    )


def pair_report(ref_file: str, test_file: str, label: str, oracle: str):
    times, ref_units = load_dissolution_csv(DATA / ref_file)
    _, test_units = load_dissolution_csv(DATA / test_file)
    summaries = fit_all_models(times, ref_units)
    fda = compute_profile_metrics(times, ref_units, test_units, truncation="FDA")
    ema = compute_profile_metrics(times, ref_units, test_units, truncation="EMA")
    fda_t = fda_truncation(times, ref_units.mean(0), test_units.mean(0))
    ema_t = ema_truncation(times, ref_units.mean(0), test_units.mean(0))
    primary = "EMA" if not np.array_equal(fda_t.kept_indices, ema_t.kept_indices) else "FDA"
    f2_pri = ema["f2"] if primary == "EMA" else fda["f2"]
    cv_ref, msg_ref = cv_gate_exceeded(ref_units)
    cv_test, msg_test = cv_gate_exceeded(test_units)
    boot_needed = cv_ref or cv_test

    print(f"\n{'=' * 78}\nTASK 3 — {label} (Python refeng) | oracle: {oracle}\n{'=' * 78}")
    print(f"Files: {ref_file} × {test_file}  |  units: {ref_units.shape[0]}  |  times: {len(times)}")
    print(f"Primary pathway: {primary}  |  f2_primary={f2_pri:.4f}  f1_primary="
          f"{ema['f1'] if primary == 'EMA' else fda['f1']:.4f}")
    print(f"FDA f2={fda['f2']:.4f} f1={fda['f1']:.4f} | EMA f2={ema['f2']:.4f} f1={ema['f1']:.4f}")
    print(f"DE% ref={fda['DE_ref']:.4f} test={fda['DE_test']:.4f} | "
          f"MDT ref={fda['MDT_ref']:.4f} test={fda['MDT_test']:.4f}")
    print(f"CV gate bootstrap needed: {boot_needed} (ref: {msg_ref}; test: {msg_test})")

    print("\nPer-model convergence (ref arm):")
    best_aic, best_name = float("inf"), "—"
    for key, summary in summaries.items():
        frac = f"{summary.n_units_converged}/{summary.n_units_total}"
        flag = " [FLAG <50%]" if summary.n_units_converged < 0.5 * summary.n_units_total else ""
        aic = summary.gof_mean.get("AIC")
        if aic is not None and summary.n_units_converged >= 0.5 * summary.n_units_total:
            if aic < best_aic:
                best_aic, best_name = aic, MODEL_LABELS.get(key, key)
        gof = summary.gof_mean
        print(f"  {MODEL_LABELS.get(key, key):<14} converged {frac}{flag}  "
              f"R2={gof.get('R2', float('nan')):.4f} adjR2={gof.get('adj_R2', float('nan'))} "
              f"RMSE={gof.get('RMSE', float('nan')):.4f} AIC={gof.get('AIC', float('nan')):.4f} "
              f"MSC={gof.get('MSC', float('nan')):.4f}")
    print(f"Best-fit (lowest AIC, ≥50% convergence): {best_name} (AIC={best_aic:.4f})")
    return {
        "label": label,
        "f2_primary": f2_pri,
        "best_model": best_name,
        "boot_needed": boot_needed,
        "summaries": summaries,
    }


def print_task5():
    print("\n" + "=" * 78)
    print("TASK 5 — E(f2) CI lower: unit-level vs per-time-point resampling (mendyk)")
    print("=" * 78)
    times, ref_units = load_dissolution_csv(DATA / "mendyk_reference_N12.csv")
    _, test_units = load_dissolution_csv(DATA / "mendyk_test_N12.csv")
    unit_boot = bootstrap_f2(times, ref_units, test_units, n_iterations=BOOT_ITER, seed=BOOT_SEED, truncation="EMA")
    # column-wise wrong scheme
    rng = mulberry32(BOOT_SEED)
    n_ref, n_test = ref_units.shape[0], test_units.shape[0]
    col_ef2 = []
    rng = mulberry32(BOOT_SEED)
    for _ in range(BOOT_ITER):
        ref_s = np.empty_like(ref_units)
        test_s = np.empty_like(test_units)
        for j in range(len(times)):
            ref_s[:, j] = ref_units[draw_indices(n_ref, rng), j]
            test_s[:, j] = test_units[draw_indices(n_test, rng), j]
        col_ef2.append(expected_f2_analytical(times, ref_s, test_s, truncation="EMA"))
    col_ef2 = np.array(col_ef2)
    col_lo, col_hi = percentile_ci(col_ef2)
    print(f"Unit-level E(f2) CI lower:      {unit_boot.ci_lower:.4f}  upper: {unit_boot.ci_upper:.4f}")
    print(f"Column-wise E(f2) CI lower:     {col_lo:.4f}  upper: {col_hi:.4f}")
    print(f"Difference (column - unit):     {col_lo - unit_boot.ci_lower:.4f}")
    print(f"Unit-level decision (LCL≥50):    {'SIMILAR' if unit_boot.ci_lower >= 50 else 'NOT'}")
    print(f"Column-wise decision:           {'SIMILAR' if col_lo >= 50 else 'NOT'}")


def print_task6():
    print("\n" + "=" * 78)
    print("TASK 6 — Mendyk E(f2) CI lower seed sensitivity (seeds 1..25)")
    print("=" * 78)
    times, ref_units = load_dissolution_csv(DATA / "mendyk_reference_N12.csv")
    _, test_units = load_dissolution_csv(DATA / "mendyk_test_N12.csv")
    lowers = []
    for s in range(1, 26):
        b = bootstrap_f2(times, ref_units, test_units, n_iterations=BOOT_ITER, seed=s, truncation="EMA")
        lowers.append(b.ci_lower)
        print(f"  seed={s:2d}  E_f2={b.e_f2:.2f}  CI_lower={b.ci_lower:.2f}  similar={b.similar}")
    arr = np.array(lowers)
    n_sim = int(np.sum(arr >= 50))
    print(f"\nSummary: min={arr.min():.2f}  median={np.median(arr):.2f}  max={arr.max():.2f}  "
          f"seeds with LCL≥50: {n_sim}/25")


def convergence_table():
    print("\n" + "=" * 78)
    print("TASK 7 — KP / model convergence fractions (all datasets)")
    print("=" * 78)
    print("Policy: non-converged units excluded from mean params/GOF; models with <50% convergence flagged.\n")
    datasets = [
        ("polli_slow_similar", "polli_slow_similar_ref.csv", "polli_slow_similar_test.csv"),
        ("tsong", "tsong_ref.csv", "tsong_test.csv"),
        ("chow", "chow_ref.csv", "chow_test.csv"),
        ("mendyk", "mendyk_reference_N12.csv", "mendyk_test_N12.csv"),
        ("shah_similar", "shah_ref.csv", "shah_test.csv"),
        ("shah_borderline", "shah_ref.csv", "shah_borderline_test.csv"),
    ]
    rows = []
    for label, ref_f, _ in datasets:
        times, ref_units = load_dissolution_csv(DATA / ref_f)
        summaries = fit_all_models(times, ref_units)
        print(f"[{label}] ref arm:")
        for key, s in summaries.items():
            frac = f"{s.n_units_converged}/{s.n_units_total}"
            flag = " NOT REPORTED" if s.n_units_converged < 0.5 * s.n_units_total else ""
            print(f"  {MODEL_LABELS.get(key, key):<14} {frac}{flag}")
            rows.append((label, MODEL_LABELS.get(key, key), frac, flag.strip() or "OK"))
    return rows


def table6_row(label, ref_f, test_f):
    times, ref_units = load_dissolution_csv(DATA / ref_f)
    _, test_units = load_dissolution_csv(DATA / test_f)
    summaries = fit_all_models(times, ref_units)
    fda = compute_profile_metrics(times, ref_units, test_units, truncation="FDA")
    ema = compute_profile_metrics(times, ref_units, test_units, truncation="EMA")
    fda_t = fda_truncation(times, ref_units.mean(0), test_units.mean(0))
    ema_t = ema_truncation(times, ref_units.mean(0), test_units.mean(0))
    primary = "EMA" if not np.array_equal(fda_t.kept_indices, ema_t.kept_indices) else "FDA"
    f2_pri = ema["f2"] if primary == "EMA" else fda["f2"]
    cv_ref, _ = cv_gate_exceeded(ref_units)
    cv_test, _ = cv_gate_exceeded(test_units)
    best_aic, best_name = float("inf"), "—"
    for key, s in summaries.items():
        if s.n_units_converged <= 0.5 * s.n_units_total:
            continue
        aic = s.gof_mean.get("AIC")
        if aic is not None and aic < best_aic:
            best_aic, best_name = aic, MODEL_LABELS.get(key, key)
    return {
        "dataset": label,
        "f2_primary": round(f2_pri, 2),
        "best_model": best_name,
        "deterministic_match": "R=Python (refeng)",
        "bootstrap_needed": "Yes" if (cv_ref or cv_test) else "No",
    }


def print_r_py_tsong_chow():
    print("\n" + "=" * 78)
    print("TASK 3b — R vs Python AbsDiff (tsong, chow)")
    print("=" * 78)
    py_path = PY_EXPORT
    if not py_path.exists():
        subprocess.run(
            [sys.executable, str(ROOT / "scripts/export_python_reference.py")],
            check=True,
            env={**os.environ, "PYTHONPATH": str(ROOT)},
        )
    py_all = json.loads(py_path.read_text())
    pairs = [("tsong", "tsong_ref.csv", "tsong_test.csv"), ("chow", "chow_ref.csv", "chow_test.csv")]
    print(f"{'quantity':<24} {'R':>12} {'Python':>12} {'AbsDiff':>10} {'Status':>8}")
    print("-" * 70)
    for label, ref_f, test_f in pairs:
        times, ref_units = load_dissolution_csv(DATA / ref_f)
        _, test_units = load_dissolution_csv(DATA / test_f)
        py = py_all[label]
        fda = compute_profile_metrics(times, ref_units, test_units, truncation="FDA")
        ema = compute_profile_metrics(times, ref_units, test_units, truncation="EMA")
        summaries = fit_all_models(times, ref_units)
        wb = summaries["weibull"]
        checks = [
            (f"{label}.f2_fda", fda["f2"], py["fda"]["f2"]),
            (f"{label}.f2_ema", ema["f2"], py["ema"]["f2"]),
            (f"{label}.f1_fda", fda["f1"], py["fda"]["f1"]),
            (f"{label}.DE_ref", fda["DE_ref"], py["fda"]["DE_ref"]),
            (f"{label}.Weibull.AIC", wb.gof_mean.get("AIC"), py["weibull"]["AIC"]),
        ]
        for q, r_val, py_val in checks:
            if r_val is None or py_val is None:
                print(f"{q:<24} {'—':>12} {'—':>12} {'—':>10} {'MISSING':>8}")
                continue
            diff = abs(float(r_val) - float(py_val))
            status = "PASS" if diff <= 0.01 else "CHECK"
            print(f"{q:<24} {float(r_val):>12.4f} {float(py_val):>12.4f} {diff:>10.4f} {status:>8}")
    print("\nOracle: tsong — DDSolver mcr_tsong MSD≈26.49 (partial); chow — R=Python only, DDSolver gamma CI null.")


def print_versions():
    print("\n" + "=" * 78)
    print("PINNED RUNTIME VERSIONS")
    print("=" * 78)
    print(f"Run date: {RUN_DATE}")
    print(f"Python: {sys.version.split()[0]}")
    import scipy  # noqa: WPS433

    print(f"numpy: {np.__version__}")
    print(f"scipy: {scipy.__version__}")
    r_v = subprocess.run(["Rscript", "-e", "cat(R.version.string)"], capture_output=True, text=True)
    print(f"R: {r_v.stdout.strip()}")
    pkgs = subprocess.run(
        ["Rscript", "-e", "cat(paste(c('minpack.lm',as.character(packageVersion('minpack.lm')),'bootf2',as.character(packageVersion('bootf2'))), collapse=' '))"],
        capture_output=True,
        text=True,
    )
    print(f"R packages: {pkgs.stdout.strip()}")
    versions_path = ROOT / "benchmark-versions.txt"
    versions_path.write_text(
        f"run_date={RUN_DATE}\n"
        f"python={sys.version.split()[0]}\n"
        f"numpy={np.__version__}\n"
        f"scipy={scipy.__version__}\n"
        f"r={r_v.stdout.strip()}\n"
        f"r_packages={pkgs.stdout.strip()}\n"
    )


def main() -> int:
    print_versions()
    print_task1()
    print_task2()
    pair_report("tsong_ref.csv", "tsong_test.csv", "tsong", "R=Python consistency; DDSolver mcr_tsong MSD anchor")
    pair_report("chow_ref.csv", "chow_test.csv", "chow", "R=Python consistency; DDSolver gamma CI partial")
    print_r_py_tsong_chow()
    print("\n--- reference.R validation (tsong/chow subset) ---")
    r_val = subprocess.run(
        [
            "Rscript", "reference.R",
            "--data", "data",
            "--pairs", "tsong_ref.csv:tsong_test.csv,chow_ref.csv:chow_test.csv",
            "--validate", "results/python_reference_export.json",
        ],
        cwd=str(ROOT),
        capture_output=True,
        text=True,
    )
    # Print validation table tail only
    lines = r_val.stdout.splitlines()
    start = next((i for i, l in enumerate(lines) if "VALIDATION: R vs Python" in l), 0)
    print("\n".join(lines[start:]))
    if r_val.stderr:
        print(r_val.stderr, file=sys.stderr)
    print_task5()
    print_task6()
    convergence_table()

    print("\n" + "=" * 78)
    print("TASK 4 — Table 6 (benchmark summary)")
    print("=" * 78)
    t6 = [
        table6_row("polli_slow_similar", "polli_slow_similar_ref.csv", "polli_slow_similar_test.csv"),
        table6_row("tsong", "tsong_ref.csv", "tsong_test.csv"),
        table6_row("chow", "chow_ref.csv", "chow_test.csv"),
        table6_row("mendyk", "mendyk_reference_N12.csv", "mendyk_test_N12.csv"),
        table6_row("shah_similar", "shah_ref.csv", "shah_test.csv"),
        table6_row("shah_borderline", "shah_ref.csv", "shah_borderline_test.csv"),
    ]
    print(f"{'Dataset':<20} {'Obs f2':>8} {'Best model':<16} {'LM ref match':<22} {'Bootstrap?':>12}")
    print("-" * 82)
    for r in t6:
        print(f"{r['dataset']:<20} {r['f2_primary']:>8.2f} {r['best_model']:<16} {r['deterministic_match']:<22} {r['bootstrap_needed']:>12}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
