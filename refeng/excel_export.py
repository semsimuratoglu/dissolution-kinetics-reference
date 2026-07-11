"""Excel workbook export for benchmark CLI results."""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import numpy as np
from openpyxl import Workbook
from openpyxl.styles import Font, PatternFill
from openpyxl.utils import get_column_letter

from refeng.bootstrap import BootstrapF2Result
from refeng.models import ModelFitSummary

HEADER_FONT = Font(bold=True)
HEADER_FILL = PatternFill(fill_type="solid", fgColor="FFD9E1F2")
NUM_FMT = "0.0000"


@dataclass
class PairExportData:
    """Structured results from one reference/test pair for reporting."""

    pair_name: str
    ref_file: str
    test_file: str
    metrics_fda: dict[str, Any]
    metrics_ema: dict[str, Any]
    ref_summaries: dict[str, ModelFitSummary] = field(default_factory=dict)
    test_summaries: dict[str, ModelFitSummary] = field(default_factory=dict)
    cv_gate_exceeded: bool = False
    cv_gate_message: str = ""
    bootstrap_required: bool = False
    bootstrap: BootstrapF2Result | None = None
    bootstrap_iterations: int = 10_000
    bootstrap_seed: int = 42


def _best_model(summaries: dict[str, ModelFitSummary]) -> tuple[str | None, float | None]:
    """Return lowest-AIC model with ≥50% unit convergence (reference arm)."""
    best_name: str | None = None
    best_aic: float | None = None
    for name, summary in summaries.items():
        if summary.n_units_converged <= 0.5 * summary.n_units_total:
            continue
        aic_val = summary.gof_mean.get("AIC")
        if aic_val is None:
            continue
        if best_aic is None or aic_val < best_aic:
            best_aic = aic_val
            best_name = name
    return best_name, best_aic


def _style_sheet(ws) -> None:
    """Apply header styling, freeze panes, and column widths."""
    if ws.max_row < 1:
        return
    for cell in ws[1]:
        cell.font = HEADER_FONT
        cell.fill = HEADER_FILL
    ws.freeze_panes = "A2"
    for col_idx in range(1, ws.max_column + 1):
        letter = get_column_letter(col_idx)
        max_len = 0
        for row in ws.iter_rows(
            min_col=col_idx, max_col=col_idx, min_row=1, max_row=ws.max_row
        ):
            val = row[0].value
            max_len = max(max_len, len(str(val)) if val is not None else 0)
        ws.column_dimensions[letter].width = min(max(max_len + 2, 10), 55)


def _write_row(ws, row_idx: int, values: list[Any], numeric_cols: set[int] | None = None) -> None:
    numeric_cols = numeric_cols or set()
    for col_idx, val in enumerate(values, start=1):
        cell = ws.cell(row=row_idx, column=col_idx, value=val)
        if col_idx in numeric_cols and isinstance(val, (int, float)):
            cell.number_format = NUM_FMT


def _summary_rows(pairs: list[PairExportData]) -> tuple[list[str], list[list[Any]], set[int]]:
    headers = [
        "Pair",
        "Reference file",
        "Test file",
        "f2 (FDA)",
        "f2 (EMA)",
        "f1 (FDA)",
        "f1 (EMA)",
        "Best-fit model (ref)",
        "Best-fit AIC (ref)",
        "DE% ref",
        "DE% test",
        "MDT ref",
        "MDT test",
        "CV gate exceeded",
        "Bootstrap needed",
        "E(f2)",
        "bc-f2",
        "90% CI lower",
        "90% CI upper",
        "Similar (y/n)",
    ]
    numeric = {4, 5, 6, 7, 9, 10, 11, 12, 13, 16, 17, 18, 19}
    rows: list[list[Any]] = []
    for pair in pairs:
        best_name, best_aic = _best_model(pair.ref_summaries)
        boot = pair.bootstrap
        rows.append(
            [
                pair.pair_name,
                pair.ref_file,
                pair.test_file,
                pair.metrics_fda["f2"],
                pair.metrics_ema["f2"],
                pair.metrics_fda["f1"],
                pair.metrics_ema["f1"],
                best_name,
                best_aic,
                pair.metrics_fda["DE_ref"],
                pair.metrics_fda["DE_test"],
                pair.metrics_fda["MDT_ref"],
                pair.metrics_fda["MDT_test"],
                "y" if pair.cv_gate_exceeded else "n",
                "y" if pair.bootstrap_required else "n",
                boot.e_f2 if boot else None,
                boot.bc_f2 if boot else None,
                boot.ci_lower if boot else None,
                boot.ci_upper if boot else None,
                ("y" if boot.similar else "n") if boot else None,
            ]
        )
    return headers, rows, numeric


ALL_PARAM_COLUMNS = ["k0", "k1", "kH", "kHC", "kKP", "n", "alpha", "beta"]
GOF_COLUMNS = ["R2", "adj_R2", "RMSE", "AIC", "MSC"]


def _param_values(params: dict[str, float]) -> list[float | None]:
    return [params.get(name) for name in ALL_PARAM_COLUMNS]


def _gof_values(gof: dict | None) -> list[float | None]:
    if not gof:
        return [None] * len(GOF_COLUMNS)
    return [gof.get(key) for key in GOF_COLUMNS]


def _models_rows(pairs: list[PairExportData]) -> tuple[list[str], list[list[Any]], set[int]]:
    headers = (
        ["Pair", "Dataset", "Model", "Unit", "Row type"]
        + ALL_PARAM_COLUMNS
        + ["R²", "adj-R²", "RMSE", "AIC", "MSC"]
    )
    param_start = 6
    gof_start = param_start + len(ALL_PARAM_COLUMNS)
    numeric = set(range(param_start, gof_start + len(GOF_COLUMNS)))
    rows: list[list[Any]] = []

    for pair in pairs:
        for dataset_label, summaries in (
            ("ref", pair.ref_summaries),
            ("test", pair.test_summaries),
        ):
            for model_key, summary in summaries.items():
                for unit_result in summary.unit_results:
                    if not unit_result.converged or not unit_result.gof:
                        continue
                    rows.append(
                        [
                            pair.pair_name,
                            dataset_label,
                            model_key,
                            unit_result.unit_index + 1,
                            "unit",
                            *_param_values(unit_result.params),
                            *_gof_values(unit_result.gof),
                        ]
                    )
                if summary.n_units_converged == 0:
                    continue
                mean_params = {
                    name: summary.param_mean[name] for name in summary.param_mean
                }
                rows.append(
                    [
                        pair.pair_name,
                        dataset_label,
                        model_key,
                        None,
                        "mean",
                        *_param_values(mean_params),
                        *_gof_values(summary.gof_mean),
                    ]
                )
                sd_params = {
                    name: summary.param_sd.get(name) for name in summary.param_mean
                }
                rows.append(
                    [
                        pair.pair_name,
                        dataset_label,
                        model_key,
                        None,
                        "sd",
                        *_param_values(sd_params),
                        *[None] * len(GOF_COLUMNS),
                    ]
                )
    return headers, rows, numeric


def _bootstrap_rows(pairs: list[PairExportData]) -> tuple[list[str], list[list[Any]], set[int]]:
    headers = [
        "Pair",
        "Iterations",
        "Seed",
        "E(f2)",
        "Median",
        "CI lower (90%)",
        "CI upper (90%)",
        "bc-f2",
        "Similar (y/n)",
        "Decision",
    ]
    numeric = {2, 4, 5, 6, 7, 8}
    rows: list[list[Any]] = []
    for pair in pairs:
        if not pair.bootstrap_required or pair.bootstrap is None:
            continue
        boot = pair.bootstrap
        decision = "similar" if boot.similar else "not similar"
        rows.append(
            [
                pair.pair_name,
                pair.bootstrap_iterations,
                pair.bootstrap_seed,
                boot.e_f2,
                boot.median,
                boot.ci_lower,
                boot.ci_upper,
                boot.bc_f2,
                "y" if boot.similar else "n",
                decision,
            ]
        )
    return headers, rows, numeric


def _golden_rows(
    golden_rows: list[dict],
) -> tuple[list[str], list[list[Any]], set[int]]:
    headers = ["Quantity", "Engine value", "Golden value", "Abs diff", "PASS/FAIL"]
    numeric = {4}
    rows: list[list[Any]] = []
    for row in golden_rows:
        rows.append(
            [
                row["quantity"],
                row["computed"],
                row["golden"],
                row["abs_diff"],
                row["status"],
            ]
        )
    return headers, rows, numeric


def write_benchmark_workbook(
    path: Path,
    pairs: list[PairExportData],
    golden_rows: list[dict] | None = None,
) -> Path:
    """Write the benchmark Excel workbook and return the resolved path."""
    path = path.resolve()
    wb = Workbook()

    # --- Summary ---
    ws_summary = wb.active
    ws_summary.title = "Summary"
    headers, rows, numeric = _summary_rows(pairs)
    _write_row(ws_summary, 1, headers)
    for i, row in enumerate(rows, start=2):
        _write_row(ws_summary, i, row, numeric)
    _style_sheet(ws_summary)

    # --- Models ---
    ws_models = wb.create_sheet("Models")
    headers, rows, numeric = _models_rows(pairs)
    _write_row(ws_models, 1, headers)
    for i, row in enumerate(rows, start=2):
        _write_row(ws_models, i, row, numeric)
    _style_sheet(ws_models)

    # --- Bootstrap ---
    ws_boot = wb.create_sheet("Bootstrap")
    headers, rows, numeric = _bootstrap_rows(pairs)
    _write_row(ws_boot, 1, headers)
    for i, row in enumerate(rows, start=2):
        _write_row(ws_boot, i, row, numeric)
    _style_sheet(ws_boot)

    # --- GoldenCheck ---
    ws_golden = wb.create_sheet("GoldenCheck")
    if golden_rows:
        headers, rows, _ = _golden_rows(golden_rows)
        numeric = {
            col
            for row in rows
            for col in (2, 3, 4)
            if col <= len(row) and isinstance(row[col - 1], (int, float))
            and not (isinstance(row[col - 1], float) and np.isnan(row[col - 1]))
        }
        _write_row(ws_golden, 1, headers)
        for i, row in enumerate(rows, start=2):
            _write_row(ws_golden, i, row, numeric)
    else:
        _write_row(ws_golden, 1, ["Quantity", "Engine value", "Golden value", "Abs diff", "PASS/FAIL"])
        _write_row(
            ws_golden,
            2,
            ["(no --golden file provided)", None, None, None, None],
        )
    _style_sheet(ws_golden)

    path.parent.mkdir(parents=True, exist_ok=True)
    wb.save(path)
    return path
