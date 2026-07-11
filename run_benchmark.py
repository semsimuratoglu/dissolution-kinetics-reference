#!/usr/bin/env python3
"""CLI entry point: fit kinetic models, compute metrics, compare to golden values."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

from refeng.bootstrap import bootstrap_f2
from refeng.excel_export import PairExportData, write_benchmark_workbook
from refeng.golden_check import compare_golden
from refeng.io import cv_gate_exceeded, load_dissolution_csv
from refeng.metrics import compute_profile_metrics
from refeng.models import fit_all_models
from refeng.polli_sets import resolve_polli_label
from refeng.shah_external_check import compare_shah_external


DETERMINISTIC_TOL = 0.01
BOOTSTRAP_TOL = 0.5

# polli_slow_similar only — do not use for fast HC demo or fast-vs-slow dissimilar pairs.
GOLDEN_PAIR_HINTS: dict[str, tuple[str, str]] = {
    "polli_slow_similar": (
        "polli_slow_similar_ref.csv",
        "polli_slow_similar_test.csv",
    ),
}


def _maybe_run_r_crosscheck(data_dir: Path, pairs_str: str) -> None:
    """Run reference.R when Rscript is available; otherwise print install hint."""
    import shutil
    import subprocess

    rscript = shutil.which("Rscript")
    if rscript is None:
        print(
            "\nR cross-check skipped (Rscript not found). "
            "Install R, then: install.packages(c('minpack.lm','bootf2'))"
        )
        return

    script = Path(__file__).resolve().parent / "reference.R"
    golden = Path(__file__).resolve().parent / "golden" / "golden_polli_slow_similar.json"
    py_export = Path(__file__).resolve().parent / "results" / "python_reference_export.json"
    cmd = [
        rscript,
        str(script),
        "--data",
        str(data_dir),
        "--pairs",
        pairs_str,
        "--validate",
        str(py_export),
        "--golden",
        str(golden),
    ]
    print("\n=== R cross-check (reference.R) ===")
    proc = subprocess.run(cmd, capture_output=False, text=True)
    if proc.returncode != 0:
        print(f"R cross-check exited with code {proc.returncode}")


def _parse_pairs(pairs_str: str) -> list[tuple[str, str]]:
    pairs: list[tuple[str, str]] = []
    for item in pairs_str.split(","):
        item = item.strip()
        if not item:
            continue
        ref, test = item.split(":")
        pairs.append((ref.strip(), test.strip()))
    return pairs


def _is_mendyk_pair(ref_file: str, test_file: str) -> bool:
    return "mendyk" in ref_file or "mendyk" in test_file


def _is_shah_pair(ref_file: str, test_file: str) -> bool:
    return ref_file == "shah_ref.csv"


def _flatten_model_results(
    pair_name: str, dataset: str, summaries: dict
) -> dict[str, float]:
    """Flatten model fit summaries into dotted keys for comparison."""
    out: dict[str, float] = {}
    for model_key, summary in summaries.items():
        prefix = f"{pair_name}.{dataset}.{model_key}"
        out[f"{prefix}.n_converged"] = float(summary.n_units_converged)
        out[f"{prefix}.n_total"] = float(summary.n_units_total)
        for pname, val in summary.param_mean.items():
            out[f"{prefix}.{pname}_mean"] = val
            out[f"{prefix}.{pname}_sd"] = summary.param_sd.get(pname, 0.0)
        for gname, val in summary.gof_mean.items():
            if val is not None:
                out[f"{prefix}.{gname}"] = val
    return out


def run_pair(
    data_dir: Path,
    ref_file: str,
    test_file: str,
    pair_name: str,
    *,
    display_name: str | None = None,
    bootstrap_iterations: int = 10_000,
    bootstrap_seed: int = 42,
) -> tuple[dict[str, float], list[str], bool, PairExportData]:
    """Run full analysis for one reference/test pair."""
    ref_path = data_dir / ref_file
    test_path = data_dir / test_file
    header = display_name or pair_name

    ref_times, ref_units = load_dissolution_csv(ref_path)
    test_times, test_units = load_dissolution_csv(test_path)

    if not np.array_equal(ref_times, test_times):
        raise ValueError(
            f"{pair_name}: reference and test time vectors must match "
            f"({ref_file} vs {test_file})"
        )

    results: dict[str, float] = {}
    log_lines: list[str] = []
    bootstrap_required = False
    ref_summaries = None
    test_summaries = None

    # --- Per-dataset kinetic fits ---
    for label, times, units in [
        ("ref", ref_times, ref_units),
        ("test", test_times, test_units),
    ]:
        summaries = fit_all_models(times, units)
        if label == "ref":
            ref_summaries = summaries
        else:
            test_summaries = summaries
        results.update(_flatten_model_results(pair_name, label, summaries))

        log_lines.append(f"\n=== {header} | {label} ({ref_file if label == 'ref' else test_file}) ===")
        for model_key, summary in summaries.items():
            log_lines.append(
                f"  {model_key}: converged {summary.n_units_converged}/{summary.n_units_total}"
                + (f" — {summary.notes}" if summary.notes else "")
            )
            if summary.param_mean:
                parts = [
                    f"{k}={summary.param_mean[k]:.6f}±{summary.param_sd.get(k, 0):.6f}"
                    for k in summary.param_mean
                ]
                log_lines.append(f"    params: {', '.join(parts)}")
            if summary.gof_mean:
                gof_parts = [
                    f"{k}={summary.gof_mean[k]:.4f}"
                    for k in summary.gof_mean
                    if summary.gof_mean[k] is not None
                ]
                log_lines.append(f"    GOF: {', '.join(gof_parts)}")

    # --- Profile metrics (FDA truncation default) ---
    metrics_fda = compute_profile_metrics(
        ref_times, ref_units, test_units, truncation="FDA"
    )
    metrics_ema = compute_profile_metrics(
        ref_times, ref_units, test_units, truncation="EMA"
    )
    for key in ("f2", "f1", "DE_ref", "DE_test", "MDT_ref", "MDT_test"):
        results[f"{pair_name}.{key}"] = metrics_fda[key]
    log_lines.append(f"\n=== {header} | similarity (FDA truncation) ===")
    log_lines.append(f"  f2={metrics_fda['f2']:.4f}, f1={metrics_fda['f1']:.4f}")
    log_lines.append(
        f"  DE% ref={metrics_fda['DE_ref']:.2f}, test={metrics_fda['DE_test']:.2f}"
    )
    log_lines.append(
        f"  MDT ref={metrics_fda['MDT_ref']:.2f}, test={metrics_fda['MDT_test']:.2f}"
    )

    # --- CV gate → bootstrap ---
    ref_cv_flag, ref_cv_msg = cv_gate_exceeded(ref_units)
    test_cv_flag, test_cv_msg = cv_gate_exceeded(test_units)
    cv_messages = []
    if ref_cv_flag:
        cv_messages.append(f"reference: {ref_cv_msg}")
    if test_cv_flag:
        cv_messages.append(f"test: {test_cv_msg}")
    cv_gate_message = "; ".join(cv_messages)

    boot = None
    if _is_shah_pair(ref_file, test_file):
        log_lines.append(f"\n=== {header} | shah: deterministic f2 anchor — bootstrap skipped ===")
    elif (ref_cv_flag or test_cv_flag) and _is_mendyk_pair(ref_file, test_file):
        bootstrap_required = True
        log_lines.append(f"\n=== {header} | CV gate EXCEEDED — bootstrap (EMA-primary) ===")
        if ref_cv_flag:
            log_lines.append(f"  reference: {ref_cv_msg}")
        if test_cv_flag:
            log_lines.append(f"  test: {test_cv_msg}")

        boot = bootstrap_f2(
            ref_times,
            ref_units,
            test_units,
            n_iterations=bootstrap_iterations,
            seed=bootstrap_seed,
            truncation="EMA",
        )
        results[f"{pair_name}.bootstrap_f2_e_f2"] = boot.e_f2
        results[f"{pair_name}.bootstrap_f2_mean"] = boot.mean
        results[f"{pair_name}.bootstrap_f2_median"] = boot.median
        results[f"{pair_name}.bootstrap_f2_ci_lower"] = boot.ci_lower
        results[f"{pair_name}.bootstrap_f2_ci_upper"] = boot.ci_upper
        results[f"{pair_name}.bootstrap_bc_f2"] = boot.bc_f2
        results[f"{pair_name}.bootstrap_similar"] = float(boot.similar)
        log_lines.append(
            f"  bootstrap: f2_obs={boot.f2_observed:.2f}, E(f2)={boot.e_f2:.2f}, "
            f"E(f2) 90% CI=[{boot.ci_lower:.2f}, {boot.ci_upper:.2f}], "
            f"bc-f2={boot.bc_f2:.2f}, similar={boot.similar} (E(f2) CI lower ≥ 50)"
        )
    else:
        log_lines.append(f"\n=== {header} | CV gate OK — bootstrap not required ===")

    export_data = PairExportData(
        pair_name=pair_name,
        ref_file=ref_file,
        test_file=test_file,
        metrics_fda=metrics_fda,
        metrics_ema=metrics_ema,
        ref_summaries=ref_summaries or {},
        test_summaries=test_summaries or {},
        cv_gate_exceeded=ref_cv_flag or test_cv_flag,
        cv_gate_message=cv_gate_message,
        bootstrap_required=bootstrap_required,
        bootstrap=boot,
        bootstrap_iterations=bootstrap_iterations,
        bootstrap_seed=bootstrap_seed,
    )

    return results, log_lines, bootstrap_required, export_data


def print_comparison_table(rows: list[dict]) -> None:
    """Print formatted comparison table."""
    header = f"{'Quantity':<55} {'Engine':>12} {'Golden':>12} {'Abs diff':>10} {'Status':>8}"
    print(header)
    print("-" * len(header))
    for row in rows:
        comp = row["computed"]
        if isinstance(comp, float):
            comp_str = f"{comp:.4f}"
        elif comp is not None:
            comp_str = str(comp)
        else:
            comp_str = "N/A"
        golden = row["golden"]
        golden_str = f"{golden:.4f}" if isinstance(golden, float) else str(golden)
        if row["abs_diff"] is not None and not (
            isinstance(row["abs_diff"], float) and np.isnan(row["abs_diff"])
        ):
            print(
                f"{row['quantity']:<55} {comp_str:>12} "
                f"{golden_str:>12} {row['abs_diff']:>10.4f} {row['status']:>8}"
            )
        else:
            print(
                f"{row['quantity']:<55} {comp_str:>12} "
                f"{golden_str:>12} {'N/A':>10} {row['status']:>8}"
            )


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Independent dissolution analytics reference benchmark"
    )
    parser.add_argument(
        "--data",
        type=str,
        default="data",
        help="Directory containing input CSV files",
    )
    parser.add_argument(
        "--pairs",
        type=str,
        required=True,
        help='Comma-separated ref:test pairs, e.g. "ref.csv:test.csv"',
    )
    parser.add_argument(
        "--golden",
        type=str,
        default=None,
        help="Path to golden expected-values JSON",
    )
    parser.add_argument(
        "--golden-shah-external",
        type=str,
        default=None,
        help="Path to external Shah literature golden JSON",
    )
    parser.add_argument(
        "--bootstrap-iterations",
        type=int,
        default=10_000,
        help="Bootstrap iteration count (default 10000)",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=42,
        help="RNG seed for bootstrap",
    )
    parser.add_argument(
        "--excel",
        type=str,
        default=None,
        help="Write results to an Excel workbook (.xlsx)",
    )
    parser.add_argument(
        "--with-r",
        action="store_true",
        help="Run reference.R cross-check when Rscript is available",
    )
    args = parser.parse_args(argv)

    data_dir = Path(args.data)
    pairs = _parse_pairs(args.pairs)

    all_results: dict[str, float] = {}
    all_logs: list[str] = []
    pair_exports: list[PairExportData] = []

    for ref_file, test_file in pairs:
        pair_name = Path(ref_file).stem + "_vs_" + Path(test_file).stem
        polli_label = resolve_polli_label(ref_file, test_file)
        display_name = polli_label or pair_name
        results, logs, _, export_data = run_pair(
            data_dir,
            ref_file,
            test_file,
            pair_name,
            display_name=display_name,
            bootstrap_iterations=args.bootstrap_iterations,
            bootstrap_seed=args.seed,
        )
        all_results.update(results)
        all_logs.extend(logs)
        pair_exports.append(export_data)

    for line in all_logs:
        print(line)

    exit_code = 0
    golden_rows: list[dict] | None = None
    if args.golden:
        golden_path = Path(args.golden)
        with golden_path.open() as fh:
            golden_data = json.load(fh)
        golden_rows, all_pass = compare_golden(pair_exports, golden_data)
        print("\n=== GoldenCheck ===")
        print_comparison_table(golden_rows)
        if not all_pass:
            exit_code = 1

    shah_external_rows: list[dict] | None = None
    if args.golden_shah_external:
        shah_path = Path(args.golden_shah_external)
        with shah_path.open() as fh:
            shah_golden = json.load(fh)
        shah_external_rows, shah_pass = compare_shah_external(pair_exports, shah_golden)
        print("\n=== ShahExternalCheck (Noce 2020 / Boddu 2024) ===")
        print_comparison_table(shah_external_rows)
        if not shah_pass:
            exit_code = 1

    if args.excel:
        excel_path = write_benchmark_workbook(
            Path(args.excel),
            pair_exports,
            golden_rows=golden_rows,
        )
        print(f"\nExcel workbook saved: {excel_path}")

    if args.with_r:
        _maybe_run_r_crosscheck(data_dir, args.pairs)

    return exit_code


if __name__ == "__main__":
    sys.exit(main())
