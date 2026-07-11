#!/usr/bin/env python3
"""Export Python reference-engine values for R cross-check validation."""

from __future__ import annotations

import json
from pathlib import Path

from refeng.bootstrap import bootstrap_f2
from refeng.io import load_dissolution_csv
from refeng.metrics import compute_profile_metrics
from refeng.models import fit_all_models

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "data"
OUT = ROOT / "results" / "python_reference_export.json"

PAIRS = {
    "polli_slow": (
        "polli_slow_similar_ref.csv",
        "polli_slow_similar_test.csv",
    ),
    "mendyk": (
        "mendyk_reference_N12.csv",
        "mendyk_test_N12.csv",
    ),
    "shah": ("shah_ref.csv", "shah_test.csv"),
    "shah_borderline": ("shah_ref.csv", "shah_borderline_test.csv"),
    "tsong": ("tsong_ref.csv", "tsong_test.csv"),
    "chow": ("chow_ref.csv", "chow_test.csv"),
}

MODEL_EXPORT = {
    "weibull": ("weibull", ["alpha", "beta"]),
    "first_order": ("fo", ["k1"]),
    "hixson_crowell": ("hc", ["kHC"]),
    "korsmeyer_peppas_full": ("kp", ["kKP", "n"]),
}


def _model_block(summaries: dict, model_key: str, export_key: str, param_names: list[str]) -> dict:
    summary = summaries[model_key]
    block: dict = {}
    for name in param_names:
        block[name] = summary.param_mean[name]
    for gof_key in ("R2", "adj_R2", "RMSE", "AIC", "MSC"):
        val = summary.gof_mean.get(gof_key)
        if val is not None:
            block[gof_key if gof_key != "adj_R2" else "adjR2"] = val
    return block


def export_pair(label: str, ref_file: str, test_file: str) -> dict:
    ref_times, ref_units = load_dissolution_csv(DATA / ref_file)
    test_times, test_units = load_dissolution_csv(DATA / test_file)
    ref_summaries = fit_all_models(ref_times, ref_units)

    fda = compute_profile_metrics(ref_times, ref_units, test_units, truncation="FDA")
    ema = compute_profile_metrics(ref_times, ref_units, test_units, truncation="EMA")

    out: dict = {
        "ref_file": ref_file,
        "test_file": test_file,
        "fda": {
            "f2": fda["f2"],
            "f1": fda["f1"],
            "DE_ref": fda["DE_ref"],
            "DE_test": fda["DE_test"],
            "MDT_ref": fda["MDT_ref"],
            "MDT_test": fda["MDT_test"],
        },
        "ema": {"f2": ema["f2"], "f1": ema["f1"]},
    }

    for model_key, (export_key, params) in MODEL_EXPORT.items():
        out[export_key] = _model_block(ref_summaries, model_key, export_key, params)

    if label == "mendyk":
        boot = bootstrap_f2(
            ref_times,
            ref_units,
            test_units,
            n_iterations=10_000,
            seed=42,
            truncation="EMA",
        )
        out["boot_ema"] = {
            "f2_observed": boot.f2_observed,
            "E_f2": boot.e_f2,
            "mean": boot.mean,
            "median": boot.median,
            "bc_f2": boot.bc_f2,
            "ci_lower": boot.ci_lower,
            "ci_upper": boot.ci_upper,
            "similar": boot.similar,
            "seed": boot.seed,
            "n_iterations": boot.n_iterations,
        }

    return out


def main() -> None:
    payload = {
        label: export_pair(label, ref_file, test_file)
        for label, (ref_file, test_file) in PAIRS.items()
    }
    OUT.write_text(json.dumps(payload, indent=2) + "\n")
    print(f"Wrote {OUT}")


if __name__ == "__main__":
    main()
