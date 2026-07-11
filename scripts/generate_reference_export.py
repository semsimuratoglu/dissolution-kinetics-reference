#!/usr/bin/env python3
"""Regenerate reference_export.md from the live reference engine."""

from __future__ import annotations

import json
import os
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

from refeng.bootstrap import bootstrap_f2
from refeng.excel_export import PairExportData
from refeng.golden_check import compare_golden
from refeng.io import cv_gate_exceeded, load_dissolution_csv
from refeng.metrics import compute_profile_metrics
from refeng.models import MODEL_SPECS, fit_model_per_unit
from refeng.polli_sets import ALL_POLLI_SETS, resolve_polli_label
from refeng.shah_external_check import compare_shah_external

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "data"
OUT_MD = ROOT / "reference_export.md"

CANONICAL_PAIRS = [
    ("polli_slow_similar_ref.csv", "polli_slow_similar_test.csv"),
    ("polli_fast_vs_slow_dissimilar_ref.csv", "polli_fast_vs_slow_dissimilar_test.csv"),
    ("mendyk_reference_N12.csv", "mendyk_test_N12.csv"),
    ("shah_ref.csv", "shah_test.csv"),
    ("shah_ref.csv", "shah_borderline_test.csv"),
]

MODEL_LABELS = {
    "zero_order": "Zero-order",
    "first_order": "First-order",
    "higuchi": "Higuchi",
    "hixson_crowell": "Hixson-Crowell",
    "korsmeyer_peppas_full": "Korsmeyer-Peppas (full profile)",
    "korsmeyer_peppas_60": "Korsmeyer-Peppas (≤60%)",
    "weibull": "Weibull",
}


def _fmt(v: float | None, digits: int = 4) -> str:
    if v is None:
        return "n/a"
    return f"{v:.{digits}f}"


def analyze_pair(ref_file: str, test_file: str) -> PairExportData:
    ref_times, ref_units = load_dissolution_csv(DATA / ref_file)
    test_times, test_units = load_dissolution_csv(DATA / test_file)
    pair_name = Path(ref_file).stem + "_vs_" + Path(test_file).stem

    ref_summaries = {
        key: fit_model_per_unit(ref_times, ref_units, key) for key in MODEL_SPECS
    }
    test_summaries = {
        key: fit_model_per_unit(test_times, test_units, key) for key in MODEL_SPECS
    }
    metrics_fda = compute_profile_metrics(ref_times, ref_units, test_units, "FDA")
    metrics_ema = compute_profile_metrics(ref_times, ref_units, test_units, "EMA")

    ref_cv, ref_msg = cv_gate_exceeded(ref_units)
    test_cv, test_msg = cv_gate_exceeded(test_units)
    boot = None
    bootstrap_required = False
    is_mendyk = "mendyk" in ref_file
    is_shah = ref_file == "shah_ref.csv"
    if is_shah:
        if ref_cv or test_cv:
            cv_note = (
                f"CV gate EXCEEDED ({ref_msg}; {test_msg}) — "
                "bootstrap skipped (Shah deterministic anchor policy)"
            )
        else:
            cv_note = "CV gate OK — bootstrap skipped (Shah deterministic anchor policy)"
    elif (ref_cv or test_cv) and is_mendyk:
        bootstrap_required = True
        boot = bootstrap_f2(
            ref_times, ref_units, test_units, n_iterations=10_000, seed=42, truncation="EMA"
        )
        cv_note = f"CV gate EXCEEDED — bootstrap EMA (ref: {ref_msg}; test: {test_msg})"
    elif ref_cv or test_cv:
        cv_note = f"CV gate exceeded but bootstrap not run ({ref_msg}; {test_msg})"
    else:
        cv_note = "CV gate OK — bootstrap not required"

    return PairExportData(
        pair_name=pair_name,
        ref_file=ref_file,
        test_file=test_file,
        metrics_fda=metrics_fda,
        metrics_ema=metrics_ema,
        ref_summaries=ref_summaries,
        test_summaries=test_summaries,
        cv_gate_exceeded=ref_cv or test_cv,
        cv_gate_message=cv_note,
        bootstrap_required=bootstrap_required,
        bootstrap=boot,
    )


def _model_table(summaries: dict, side: str) -> list[str]:
    lines = [
        f"#### {side} profile — model fits",
        "",
        "| Model | Convergence | Parameters (mean±SD) | R² | adj-R² | RMSE | AIC | MSC |",
        "|-------|-------------|----------------------|-----|--------|------|-----|-----|",
    ]
    for key, label in MODEL_LABELS.items():
        sm = summaries[key]
        if sm.n_units_converged == 0:
            params = "n/a"
            gof = ("n/a",) * 5
        else:
            parts = [
                f"{p}={sm.param_mean[p]:.6f}±{sm.param_sd.get(p, 0):.6f}"
                for p in sm.param_mean
            ]
            params = ", ".join(parts)
            gof = (
                _fmt(sm.gof_mean.get("R2")),
                _fmt(sm.gof_mean.get("adj_R2")),
                _fmt(sm.gof_mean.get("RMSE")),
                _fmt(sm.gof_mean.get("AIC")),
                _fmt(sm.gof_mean.get("MSC")),
            )
        lines.append(
            f"| {label} | {sm.n_units_converged}/{sm.n_units_total} | {params} | "
            f"{gof[0]} | {gof[1]} | {gof[2]} | {gof[3]} | {gof[4]} |"
        )
        if sm.notes:
            lines.append(f"| | *{sm.notes}* | | | | | | |")
    lines.append("")
    return lines


def _pair_section(pair: PairExportData) -> list[str]:
    label = resolve_polli_label(pair.ref_file, pair.test_file) or pair.pair_name
    ref_times, _ = load_dissolution_csv(DATA / pair.ref_file)
    _, ref_units = load_dissolution_csv(DATA / pair.ref_file)
    _, test_units = load_dissolution_csv(DATA / pair.test_file)
    lines = [
        f"### {label} (`{pair.ref_file}` vs `{pair.test_file}`)",
        "",
        f"- Time points: {len(ref_times)} | Ref units: {ref_units.shape[0]} | "
        f"Test units: {test_units.shape[0]}",
        f"- Time vector: {', '.join(f'{t:g}' for t in ref_times)}",
        "",
    ]
    lines.extend(_model_table(pair.ref_summaries, "Reference"))
    lines.extend(_model_table(pair.test_summaries, "Test"))
    m_fda, m_ema = pair.metrics_fda, pair.metrics_ema
    lines.extend(
        [
            "#### Similarity (reference mean vs test mean)",
            "",
            "| Metric | FDA | EMA |",
            "|--------|-----|-----|",
            f"| f₂ | {_fmt(m_fda['f2'])} | {_fmt(m_ema['f2'])} |",
            f"| f₁ | {_fmt(m_fda['f1'])} | {_fmt(m_ema['f1'])} |",
            f"| DE% (ref) | {_fmt(m_fda['DE_ref'], 2)} | n/a |",
            f"| DE% (test) | {_fmt(m_fda['DE_test'], 2)} | n/a |",
            f"| MDT (ref) | {_fmt(m_fda['MDT_ref'], 2)} | n/a |",
            f"| MDT (test) | {_fmt(m_fda['MDT_test'], 2)} | n/a |",
            f"| Truncation times | {', '.join(f'{t:g}' for t in m_fda['truncated_times'])} | "
            f"{', '.join(f'{t:g}' for t in m_ema['truncated_times'])} |",
            "",
            "#### CV gate / Bootstrap",
            "",
            f"**{pair.cv_gate_message}**",
            "",
        ]
    )
    if pair.bootstrap is not None:
        b = pair.bootstrap
        lines.extend(
            [
                "| Field | Value |",
                "|-------|-------|",
                f"| f₂ observed (EMA) | {_fmt(b.f2_observed)} |",
                f"| E(f₂) point | {_fmt(b.e_f2)} |",
                f"| E(f₂) median (bootstrap) | {_fmt(b.median)} |",
                f"| bc-f₂ | {_fmt(b.bc_f2)} |",
                f"| 90% CI lower | {_fmt(b.ci_lower)} |",
                f"| 90% CI upper | {_fmt(b.ci_upper)} |",
                f"| similar (CI lower ≥ 50) | {b.similar} |",
                f"| seed | {b.seed} |",
                "",
            ]
        )
    return lines


def main() -> None:
    pairs = [analyze_pair(ref, test) for ref, test in CANONICAL_PAIRS]
    now = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC")

    polli_golden = json.loads((ROOT / "golden/golden_polli_slow_similar.json").read_text())
    shah_golden = json.loads((ROOT / "golden/golden_shah_external.json").read_text())
    polli_rows, polli_pass = compare_golden(
        [p for p in pairs if resolve_polli_label(p.ref_file, p.test_file) == "polli_slow_similar"],
        polli_golden,
    )
    shah_rows, shah_pass = compare_shah_external(
        [p for p in pairs if p.ref_file == "shah_ref.csv"],
        shah_golden,
    )

    lines = [
        "# Reference engine export",
        "",
        "Full export from the independent Python reference engine in this repository.",
        f"Generated: {now}",
        "",
        "## Canonical dataset registry",
        "",
        "| Label | Ref | Test | Role |",
        "|-------|-----|------|------|",
    ]
    for ps in ALL_POLLI_SETS:
        test = ps.test_file or "—"
        lines.append(f"| `{ps.label}` | `{ps.ref_file}` | `{test}` | {ps.role} |")
    lines.extend(
        [
            "| `mendyk_high_cv` | `mendyk_reference_N12.csv` | `mendyk_test_N12.csv` | Synthetic high-CV bootstrap anchor |",
            "| `shah_similar` | `shah_ref.csv` | `shah_test.csv` | bootf2::shah1998 test1, f₂≈60.03 |",
            "| `shah_borderline` | `shah_ref.csv` | `shah_borderline_test.csv` | bootf2::shah1998 test5, f₂≈48.05 |",
            "",
            "---",
            "",
            "## 1. golden_polli_slow_similar.json",
            "",
            "```json",
            json.dumps(polli_golden, indent=2),
            "```",
            "",
            "## 2. golden_shah_external.json (Noce 2020 / Boddu 2024)",
            "",
            "```json",
            json.dumps(shah_golden, indent=2),
            "```",
            "",
            "## 3. Python reference engine — canonical pairs",
            "",
            "Conventions: per-unit fit → mean±SD; RMSE = √(RSS/(n−p)); AIC = N·ln(RSS)+2P;",
            "f₂/f₁ with FDA+EMA 85% truncation; DE%/MDT on full profile. Bootstrap for mendyk only (EMA-primary).",
            "",
        ]
    )
    for pair in pairs:
        lines.extend(_pair_section(pair))

    lines.extend(
        [
            "---",
            "",
            "## 4. Golden / external validation",
            "",
            f"**polli_slow_similar GoldenCheck:** {'PASS' if polli_pass else 'CHECK'}",
            "",
            "| Quantity | Engine | Golden | Abs diff | Status |",
            "|----------|--------|--------|----------|--------|",
        ]
    )
    for row in polli_rows:
        comp = row["computed"]
        comp_s = f"{comp:.4f}" if isinstance(comp, float) else str(comp)
        golden_s = f"{row['golden']:.4f}" if isinstance(row['golden'], float) else str(row["golden"])
        diff_s = f"{row['abs_diff']:.4g}" if row["abs_diff"] is not None else "N/A"
        lines.append(
            f"| {row['quantity']} | {comp_s} | {golden_s} | {diff_s} | {row['status']} |"
        )

    lines.extend(
        [
            "",
            f"**ShahExternalCheck (literature):** {'PASS' if shah_pass else 'CHECK'}",
            "",
            "| Quantity | Engine | Golden | Abs diff | Status |",
            "|----------|--------|--------|----------|--------|",
        ]
    )
    for row in shah_rows:
        comp = row["computed"]
        comp_s = f"{comp:.4f}" if isinstance(comp, float) else str(comp)
        golden_s = f"{row['golden']:.4f}" if isinstance(row['golden'], float) else str(row["golden"])
        diff_s = f"{row['abs_diff']:.4g}" if row["abs_diff"] is not None else "N/A"
        lines.append(
            f"| {row['quantity']} | {comp_s} | {golden_s} | {diff_s} | {row['status']} |"
        )

    lines.extend(
        [
            "",
            "## 5. Shah borderline summary",
            "",
        ]
    )
    borderline = next(p for p in pairs if p.test_file == "shah_borderline_test.csv")
    lines.extend(
        [
            f"- **f₂ (FDA/EMA):** {_fmt(borderline.metrics_fda['f2'])} (literature ≈48.05)",
            f"- **CV gate:** {borderline.cv_gate_message}",
            f"- **Bootstrap triggered:** {borderline.bootstrap_required}",
            "- **Verdict:** borderline similarity (point estimate near 50; literature percentile CI upper < 50 for test5)",
            "",
            "## 6. R cross-check",
            "",
            "Run:",
            "",
            "```bash",
            "Rscript reference.R --data data \\",
            '  --pairs "polli_slow_similar_ref.csv:polli_slow_similar_test.csv,'
            "mendyk_reference_N12.csv:mendyk_test_N12.csv,shah_ref.csv:shah_test.csv,shah_ref.csv:shah_borderline_test.csv\" \\",
            "  --validate results/python_reference_export.json \\",
            "  --golden golden/golden_polli_slow_similar.json",
            "```",
            "",
            "Shah export: `Rscript scripts/export_shah1998.R --out-dir data`",
            "",
        ]
    )

    OUT_MD.write_text("\n".join(lines) + "\n")
    print(f"Wrote {OUT_MD} ({len(lines)} lines)")

    subprocess.run(
        [sys.executable, str(ROOT / "scripts/export_python_reference.py")],
        check=True,
        env={**os.environ, "PYTHONPATH": str(ROOT)},
    )


if __name__ == "__main__":
    main()
