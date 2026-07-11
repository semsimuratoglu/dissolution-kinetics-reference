"""Compare Shah pairs against externally published literature anchors."""

from __future__ import annotations

from typing import Any

from refeng.excel_export import PairExportData


def _compare_value(
    quantity: str,
    computed: float | None,
    golden_val: float,
    tol: float,
) -> dict[str, Any]:
    if computed is None:
        return {
            "quantity": quantity,
            "computed": None,
            "golden": golden_val,
            "abs_diff": None,
            "status": "MISSING",
        }
    abs_diff = abs(computed - golden_val)
    status = "PASS" if abs_diff <= tol else "CHECK"
    return {
        "quantity": quantity,
        "computed": computed,
        "golden": golden_val,
        "abs_diff": abs_diff,
        "status": status,
    }


def _resolve_shah_section(pair: PairExportData, golden: dict[str, Any]) -> str | None:
    ref = pair.ref_file
    test = pair.test_file
    if ref == "shah_ref.csv" and test == "shah_test.csv":
        return "shah_similar"
    if ref == "shah_ref.csv" and test == "shah_borderline_test.csv":
        return "shah_borderline"
    for key, section in golden.items():
        if key.startswith("_"):
            continue
        if section.get("ref_file") == ref and section.get("test_file") == test:
            return key
    return None


def compare_shah_external_pair(
    pair: PairExportData,
    golden: dict[str, Any],
) -> list[dict[str, Any]]:
    """Compare a Shah pair against golden_shah_external.json."""
    section_key = _resolve_shah_section(pair, golden)
    if section_key is None:
        return []

    meta = golden.get("_meta", {})
    tol_f2 = float(meta.get("tolerance_f2", 0.01))
    tol_ci = float(meta.get("tolerance_ci", 2.0))
    section = golden[section_key]
    prefix = section.get("pair_label", section_key)
    rows: list[dict[str, Any]] = []

    metrics = pair.metrics_fda
    rows.append(
        _compare_value(
            f"{prefix}.f2_observed",
            float(metrics["f2"]),
            float(section["f2_observed"]),
            tol_f2,
        )
    )
    rows.append(
        _compare_value(
            f"{prefix}.f1_observed",
            float(metrics["f1"]),
            float(section["f1_observed"]),
            tol_f2,
        )
    )

    boot = section.get("bootstrap", {})
    if pair.bootstrap is not None:
        rows.append(
            _compare_value(
                f"{prefix}.boot.percentile_ci_lower",
                float(pair.bootstrap.ci_lower),
                float(boot["percentile_90_ci_lower"]),
                tol_ci,
            )
        )
        rows.append(
            _compare_value(
                f"{prefix}.boot.percentile_ci_upper",
                float(pair.bootstrap.ci_upper),
                float(boot["percentile_90_ci_upper"]),
                tol_ci,
            )
        )
        rows.append(
            _compare_value(
                f"{prefix}.boot.bc_f2",
                float(pair.bootstrap.bc_f2),
                float(section["f2_observed"]),
                tol_ci,
            )
        )
    else:
        for key in (
            "percentile_90_ci_lower",
            "percentile_90_ci_upper",
        ):
            if key in boot:
                rows.append(
                    {
                        "quantity": f"{prefix}.boot.{key}",
                        "computed": None,
                        "golden": float(boot[key]),
                        "abs_diff": None,
                        "status": "SKIP_NO_BOOT",
                    }
                )

    expected_similar = boot.get("similar_by_percentile")
    if pair.bootstrap is not None and expected_similar is not None:
        computed_similar = pair.bootstrap.similar
        rows.append(
            {
                "quantity": f"{prefix}.boot.similar",
                "computed": str(computed_similar),
                "golden": str(expected_similar),
                "abs_diff": None,
                "status": "PASS" if computed_similar == expected_similar else "CHECK",
            }
        )

    return rows


def compare_shah_external(
    pairs: list[PairExportData],
    golden: dict[str, Any],
) -> tuple[list[dict[str, Any]], bool]:
    """Run external Shah literature checks for matching pairs."""
    rows: list[dict[str, Any]] = []
    for pair in pairs:
        rows.extend(compare_shah_external_pair(pair, golden))
    actionable = [r for r in rows if r["status"] not in ("SKIP_NO_BOOT",)]
    all_pass = all(r["status"] == "PASS" for r in actionable)
    return rows, all_pass
