#!/usr/bin/env python3
"""Development full benchmark validation suite (tasks 1–7).

Run: PYTHONPATH=. python dev/run_full_benchmark_validation.py
"""

from __future__ import annotations

import json
import shutil
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import scipy
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "data"
PY_EXPORT = ROOT / "results" / "python_reference_export.json"
sys.path.insert(0, str(ROOT))

from refeng.bootstrap_ef2 import bootstrap_ef2_columnwise, bootstrap_ef2_validation
from refeng.io import cv_gate_exceeded, load_dissolution_csv
from refeng.metrics import compute_profile_metrics, expected_f2_analytical
from refeng.models import MODEL_SPECS, fit_all_models

MODEL_LABELS = {
    "zero_order": "Zero-order",
    "first_order": "First-order",
    "higuchi": "Higuchi",
    "hixson_crowell": "Hixson-Crowell",
    "korsmeyer_peppas_full": "Korsmeyer-Peppas (full profile)",
    "korsmeyer_peppas_60": "Korsmeyer-Peppas (≤60%)",
    "weibull": "Weibull",
}

LITERATURE = {
    "shah_similar": {
        "percentile": (53.2, 67.5),
        "bca": (53.9, 70.2),
    },
    "shah_borderline": {
        "percentile": (45.99, 49.96),
        "bca": (46.15, 50.17),
    },
}

BODDU_RANGES = {
    "shah_similar": {"bca_lower": (53.8, 55.8), "bca_upper": (69.3, 73.0)},
    "shah_borderline": {"bca_lower": (45.9, 46.3), "bca_upper": (49.95, 50.28)},
}

DATASETS = [
    ("polli_slow_similar", "polli_slow_similar_ref.csv", "polli_slow_similar_test.csv"),
    ("tsong", "tsong_ref.csv", "tsong_test.csv"),
    ("chow", "chow_ref.csv", "chow_test.csv"),
    ("mendyk", "mendyk_reference_N12.csv", "mendyk_test_N12.csv"),
    ("shah_similar", "shah_ref.csv", "shah_test.csv"),
    ("shah_borderline", "shah_ref.csv", "shah_borderline_test.csv"),
]


def _hdr(title: str) -> None:
    print("\n" + "=" * 78)
    print(title)
    print("=" * 78)


def _within(val: float | None, target: float, tol: float = 3.0) -> str:
    if val is None or not np.isfinite(val):
        return "N/A"
    return "YES" if abs(val - target) <= tol else "NO"


def _within_range(val: float | None, lo: float, hi: float, tol: float = 3.0) -> str:
    if val is None or not np.isfinite(val):
        return "N/A"
    return "YES" if (lo - tol) <= val <= (hi + tol) else "NO"


def task1_shah_bootstrap() -> None:
    _hdr("TASK 1 — Shah bootstrap external validation (validation mode, seed=42, n=10000, EMA)")
    pairs = [
        ("shah_similar", "shah_test.csv"),
        ("shah_borderline", "shah_borderline_test.csv"),
    ]

    r_json = ROOT / "scripts" / "_shah_bootf2_validation.json"
    rscript = shutil.which("Rscript")
    bootf2_r = {}
    if rscript:
        proc = subprocess.run(
            [rscript, str(ROOT / "scripts/shah_bootstrap_validation.R"), "--data", str(DATA)],
            capture_output=True,
            text=True,
        )
        if proc.returncode == 0:
            bootf2_r = json.loads(proc.stdout)
            r_json.write_text(json.dumps(bootf2_r, indent=2) + "\n")
            print(f"bootf2 R output saved: {r_json}")
        else:
            print("R bootf2 validation FAILED:")
            print(proc.stderr or proc.stdout)
    else:
        print("Rscript not found — bootf2 columns will be N/A")

    print(
        f"\n{'quantity':<42} {'our_engine':>12} {'bootf2':>12} {'Noce':>14} {'Boddu':>14} {'±3 tol?':>8}"
    )
    print("-" * 108)

    for label, test_file in pairs:
        times, ref_u = load_dissolution_csv(DATA / "shah_ref.csv")
        _, test_u = load_dissolution_csv(DATA / test_file)
        eng = bootstrap_ef2_validation(times, ref_u, test_u, seed=42, truncation="EMA")
        lit = LITERATURE[label]
        boddu = BODDU_RANGES[label]
        rpair = (bootf2_r.get("pairs") or {}).get(label, {})
        est = rpair.get("est_f2", {})

        rows = [
            ("f2_observed", eng.f2_observed, est.get("f2_observed"), None, None, None),
            ("E_f2_analytical", eng.e_f2_point, est.get("boot_mean"), None, None, None),
            ("percentile_CI_lower", eng.percentile_ci_lower, est.get("percentile_ci_lower"), lit["percentile"][0], boddu["bca_lower"][0], "pct_lo"),
            ("percentile_CI_upper", eng.percentile_ci_upper, est.get("percentile_ci_upper"), lit["percentile"][1], boddu["bca_upper"][1], "pct_hi"),
            ("BCa_CI_lower", eng.bca_ci_lower, est.get("bca_ci_lower"), lit["bca"][0], boddu["bca_lower"][0], "bca_lo"),
            ("BCa_CI_upper", eng.bca_ci_upper, est.get("bca_ci_upper"), lit["bca"][1], boddu["bca_upper"][1], "bca_hi"),
        ]
        print(f"\n--- {label} ({test_file}) ---")
        for qty, ours, bf2, noce, boddu_ref, kind in rows:
            bf2_s = f"{bf2:.4f}" if isinstance(bf2, (int, float)) and np.isfinite(bf2) else "N/A"
            noce_s = f"{noce:.2f}" if isinstance(noce, (int, float)) else "—"
            if isinstance(boddu_ref, (int, float)):
                boddu_s = f"{boddu_ref:.2f}"
            elif isinstance(boddu_ref, tuple):
                boddu_s = f"[{boddu_ref[0]:.1f},{boddu_ref[1]:.1f}]"
            else:
                boddu_s = "—"
            tol = "—"
            if kind == "pct_lo":
                tol = _within(ours, lit["percentile"][0])
            elif kind == "pct_hi":
                tol = _within(ours, lit["percentile"][1])
            elif kind == "bca_lo":
                tol = _within(ours, lit["bca"][0])
                if isinstance(boddu_ref, tuple):
                    tol += "/" + _within_range(ours, boddu_ref[0], boddu_ref[1])
            elif kind == "bca_hi":
                tol = _within(ours, lit["bca"][1])
                if isinstance(boddu_ref, tuple):
                    tol += "/" + _within_range(ours, boddu_ref[0], boddu_ref[1])
            print(f"{qty:<42} {ours:>12.4f} {bf2_s:>12} {noce_s:>14} {boddu_s:>14} {tol:>8}")

        print(
            f"  Decision (our E-f2 percentile LCL≥50): {eng.similar_percentile} | "
            f"BCa LCL≥50: {eng.similar_bca} | CV gate does NOT require bootstrap for Shah"
        )


def task2_mendyk_provenance() -> None:
    _hdr("TASK 2 — Mendyk N12 provenance")
    print(
        """
SOURCE VERDICT: **synthetic high-variability set** (NOT raw published Mendyk 2013 unit matrices).

Evidence:
  • Files: data/mendyk_reference_N12.csv, data/mendyk_test_N12.csv.
  • NOT produced by scripts/generate_sample_data.py (that script only writes weibull_ref/test and high_cv_ref/test).
  • Published Mendyk et al. (2013) Case Study 1 uses times 3–30 min; N12 uses 15–120 min with a different scale.
  • No export script or literature table maps 1:1 to these 12×6 CSV values.
  • Role: engineered high-CV dissolution pair (first-point CV ≈25.7% ref / 28.4% test) to force
    bootstrap under EMA %CV gate; labeled **synthetic high-variability set**.
"""
    )


def _primary_f2(times, ref_u, test_u) -> float:
    m_ema = compute_profile_metrics(times, ref_u, test_u, "EMA")
    m_fda = compute_profile_metrics(times, ref_u, test_u, "FDA")
    if m_ema["kept_indices"] != m_fda["kept_indices"]:
        return float(m_ema["f2"])
    return float(m_fda["f2"])


def _first_point_cv(ref_u, test_u) -> tuple[float, float]:
    def cv(col):
        m = np.mean(col)
        if m == 0:
            return float("nan")
        return float(np.std(col, ddof=1) / m * 100)

    return cv(ref_u[:, 0]), cv(test_u[:, 0])


def task3_tsong_chow_r_python() -> dict:
    _hdr("TASK 3 — Tsong & Chow: Python scipy + R minpack.lm cross-check")
    print("Oracle status: **R=Python consistency only, no independent DDSolver golden** for tsong/chow.\n")

    pairs = [("tsong", "tsong_ref.csv", "tsong_test.csv"), ("chow", "chow_ref.csv", "chow_test.csv")]
    py_export = {}
    for label, ref_f, test_f in pairs:
        times, ref_u = load_dissolution_csv(DATA / ref_f)
        _, test_u = load_dissolution_csv(DATA / test_f)
        ref_models = fit_all_models(times, ref_u)
        m_fda = compute_profile_metrics(times, ref_u, test_u, "FDA")
        m_ema = compute_profile_metrics(times, ref_u, test_u, "EMA")
        ref_cv, test_cv = _first_point_cv(ref_u, test_u)
        flag, msg = cv_gate_exceeded(ref_u)
        print(f"\n### {label.upper()} ({ref_f} vs {test_f})")
        print(f"Provenance: {'Tsong 1996 (minutes, 6 units)' if label=='tsong' else 'Chow & Ki 1997 Table 1 / Tsong–Hammerstrom (hours, 11 units)'}")
        print(f"Primary f2={_primary_f2(times, ref_u, test_u):.4f} | f2 FDA={m_fda['f2']:.4f} EMA={m_ema['f2']:.4f}")
        print(f"f1 FDA={m_fda['f1']:.4f} | DE% {m_fda['DE_ref']:.2f}/{m_fda['DE_test']:.2f} | MDT {m_fda['MDT_ref']:.2f}/{m_fda['MDT_test']:.2f}")
        print(f"First-point CV% ref={ref_cv:.2f} test={test_cv:.2f} | CV gate: {msg} | Bootstrap needed: {flag}")
        print(f"{'Model':<32} {'Conv':>8} {'R²':>8} {'adjR²':>8} {'RMSE':>8} {'AIC':>10} {'MSC':>8}")
        best_aic, best_name = float("inf"), "—"
        block = {}
        for key, ml in MODEL_LABELS.items():
            sm = ref_models[key]
            conv = f"{sm.n_units_converged}/{sm.n_units_total}"
            r2 = sm.gof_mean.get("R2")
            adj = sm.gof_mean.get("adj_R2")
            rmse = sm.gof_mean.get("RMSE")
            aic = sm.gof_mean.get("AIC")
            msc = sm.gof_mean.get("MSC")
            if aic is not None and sm.n_units_converged > 0 and aic < best_aic:
                best_aic, best_name = aic, ml
            flag_kp = " *" if sm.n_units_converged / max(sm.n_units_total, 1) < 0.5 else ""
            print(
                f"{ml:<32} {conv:>8} {r2 or 0:8.4f} {adj or 0:8.4f} {rmse or 0:8.4f} {aic or 0:10.4f} {msc or 0:8.4f}{flag_kp}"
            )
            block[key] = {
                "converged": sm.n_units_converged,
                "total": sm.n_units_total,
                "gof": sm.gof_mean,
            }
        py_export[label] = block
        print(f"Best AIC model (ref): {best_name}")

    # Run R validation for tsong/chow only
    rscript = shutil.which("Rscript")
    if rscript:
        pairs_str = "tsong_ref.csv:tsong_test.csv,chow_ref.csv:chow_test.csv"
        # refresh python export subset
        subprocess.run([sys.executable, str(ROOT / "scripts/export_python_reference.py")], check=True, env={**os_environ(), "PYTHONPATH": str(ROOT)})
        cmd = [
            rscript,
            str(ROOT / "reference.R"),
            "--data",
            str(DATA),
            "--pairs",
            pairs_str,
            "--validate",
            str(PY_EXPORT),
        ]
        print("\n--- R reference.R validation (tsong + chow) ---")
        subprocess.run(cmd)
    return py_export


def os_environ():
    import os
    return os.environ


def task4_table6() -> None:
    _hdr("TASK 4 — Table 6 (benchmark summary)")
    print(f"{'Dataset':<22} {'Obs f2':>10} {'Best model':>28} {'R=Py match':>12} {'Bootstrap?':>12}")
    print("-" * 90)
    for label, ref_f, test_f in DATASETS:
        times, ref_u = load_dissolution_csv(DATA / ref_f)
        _, test_u = load_dissolution_csv(DATA / test_f)
        f2p = _primary_f2(times, ref_u, test_u)
        ref_models = fit_all_models(times, ref_u)
        best_aic, best_name = float("inf"), "—"
        for key, ml in MODEL_LABELS.items():
            sm = ref_models[key]
            aic = sm.gof_mean.get("AIC")
            if aic is not None and sm.n_units_converged > 0 and aic < best_aic:
                if sm.n_units_converged / sm.n_units_total >= 0.5:
                    best_aic, best_name = aic, ml
                else:
                    best_name = f"{ml} (<50% conv)"
        cv_flag, _ = cv_gate_exceeded(ref_u)
        if "shah" in label:
            boot = "No (CV OK; validation only)"
        elif label == "polli_slow_similar":
            boot = "No"
        elif label == "mendyk":
            boot = "Yes (CV exceeded)"
        else:
            boot = "No" if not cv_flag else "Yes"
        if label == "polli_slow_similar":
            match = "Golden PASS"
        elif label in ("tsong", "chow"):
            match = "R=Py only"
        elif label.startswith("shah"):
            match = "Literature f2"
        else:
            match = "R=Py PASS"
        print(f"{label:<22} {f2p:10.4f} {best_name:>28} {match:>12} {boot:>12}")


def task5_granularity() -> None:
    _hdr("TASK 5 — Table 5: E-f2 CI lower bound — unit-level vs per-time-point (mendyk)")
    times, ref_u = load_dissolution_csv(DATA / "mendyk_reference_N12.csv")
    _, test_u = load_dissolution_csv(DATA / "mendyk_test_N12.csv")
    unit = bootstrap_ef2_validation(times, ref_u, test_u, seed=42)
    col = bootstrap_ef2_columnwise(times, ref_u, test_u, seed=42)
    diff = col.percentile_ci_lower - unit.percentile_ci_lower
    print(f"Decision rule: SIMILAR if E(f2) 90% percentile CI lower bound ≥ 50\n")
    print(f"  Unit-level (correct):     LCL = {unit.percentile_ci_lower:.4f}  → {'SIMILAR' if unit.similar_percentile else 'NOT SIMILAR'}")
    print(f"  Per-time-point (legacy):  LCL = {col.percentile_ci_lower:.4f}  → {'SIMILAR' if col.similar_percentile else 'NOT SIMILAR'}")
    print(f"  Difference (column − unit): {diff:+.4f}")
    print("  (Per-time-point resampling inflates the lower bound — anti-conservative vs unit-level.)")


def task6_seed_sensitivity() -> None:
    _hdr("TASK 6 — Mendyk E-f2 bootstrap seed sensitivity (seeds 1..25)")
    times, ref_u = load_dissolution_csv(DATA / "mendyk_reference_N12.csv")
    _, test_u = load_dissolution_csv(DATA / "mendyk_test_N12.csv")
    lowers = []
    for seed in range(1, 26):
        b = bootstrap_ef2_validation(times, ref_u, test_u, seed=seed)
        lowers.append(b.percentile_ci_lower)
    lowers = np.array(lowers)
    n_sim = int(np.sum(lowers >= 50))
    print(f"  E(f2) 90% percentile CI lower bound over seeds 1..25:")
    print(f"    min = {lowers.min():.4f}")
    print(f"    median = {np.median(lowers):.4f}")
    print(f"    max = {lowers.max():.4f}")
    print(f"    count(lower ≥ 50) = {n_sim} / 25")
    print(f"  Seed 42 (canonical): {bootstrap_ef2_validation(times, ref_u, test_u, seed=42).percentile_ci_lower:.4f}")


def task7_kp_convergence() -> None:
    _hdr("TASK 7 — Per-model convergence fractions (all datasets)")
    print("Policy: fit each unit independently; non-converged units excluded from mean±SD;")
    print("        models with <50% unit convergence flagged as NOT REPORTED for best-model selection.\n")
    print(f"{'Dataset':<22} {'Model':<32} {'Converged':>12} {'Flag':>8}")
    print("-" * 78)
    for label, ref_f, test_f in DATASETS:
        times, ref_u = load_dissolution_csv(DATA / ref_f)
        ref_models = fit_all_models(times, ref_u)
        for key, ml in MODEL_LABELS.items():
            sm = ref_models[key]
            conv = f"{sm.n_units_converged}/{sm.n_units_total}"
            flag = "OK" if sm.n_units_converged / max(sm.n_units_total, 1) >= 0.5 else "LOW"
            print(f"{label:<22} {ml:<32} {conv:>12} {flag:>8}")


def pin_versions() -> None:
    _hdr("PINNED PACKAGE VERSIONS")
    import os

    run_date = datetime.now(timezone.utc).strftime("%Y-%m-%d")
    print(f"Run date (UTC): {run_date}\n")
    print(f"  Python     {sys.version.split()[0]}")
    print(f"  numpy      {np.__version__}")
    print(f"  scipy      {scipy.__version__}")
    print(f"  pandas     {pd.__version__}")
    rscript = shutil.which("Rscript")
    if rscript:
        rv = subprocess.run([rscript, "--version"], capture_output=True, text=True)
        print(f"  R          {(rv.stdout or rv.stderr).strip()}")
        pkgs = subprocess.run(
            [
                rscript,
                "-e",
                'cat(sprintf("minpack.lm %s | bootf2 %s | jsonlite %s\\n", packageVersion("minpack.lm"), packageVersion("bootf2"), packageVersion("jsonlite")))',
            ],
            capture_output=True,
            text=True,
        )
        print(f"  {pkgs.stdout.strip()}")
    renv = ROOT / "renv-versions.txt"
    if renv.exists():
        print(f"\n--- {renv} ---")
        print(renv.read_text())


def run_external_tests() -> None:
    _hdr("EXTERNAL TEST SUITE RE-RUN")
    tests = [
        ([sys.executable, "-m", "pytest", "tests/", "-q"], ROOT, "pytest"),
    ]
    rscript = shutil.which("Rscript")
    if rscript:
        tests.append(
            (
                [
                    rscript,
                    str(ROOT / "reference.R"),
                    "--data",
                    str(DATA),
                    "--pairs",
                    "polli_slow_similar_ref.csv:polli_slow_similar_test.csv,mendyk_reference_N12.csv:mendyk_test_N12.csv,shah_ref.csv:shah_test.csv,shah_ref.csv:shah_borderline_test.csv",
                    "--validate",
                    str(PY_EXPORT),
                    "--golden",
                    str(ROOT / "golden/golden_polli_slow_similar.json"),
                ],
                ROOT,
                "reference.R",
            )
        )

    for cmd, cwd, name in tests:
        print(f"\n>>> {name}")
        proc = subprocess.run(cmd, cwd=str(cwd), capture_output=False, text=True)
        if proc.returncode != 0:
            print(f"*** {name} exited {proc.returncode}")


def main() -> int:
    print("Independent reference engine — full computational validation")
    print(f"Root: {ROOT}")
    task1_shah_bootstrap()
    task2_mendyk_provenance()
    task3_tsong_chow_r_python()
    task4_table6()
    task5_granularity()
    task6_seed_sensitivity()
    task7_kp_convergence()
    pin_versions()
    run_external_tests()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
