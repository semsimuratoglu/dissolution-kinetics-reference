"""Golden-value comparison for externally anchored benchmark JSON."""

from __future__ import annotations

from typing import Any

from refeng.excel_export import PairExportData
from refeng.models import ModelFitSummary
from refeng.polli_sets import POLLI_SLOW_SIMILAR, golden_section_key, is_golden_polli_pair

# Golden display name -> engine model key
GOLDEN_MODEL_KEYS: dict[str, str] = {
    "Weibull": "weibull",
    "First-order": "first_order",
    "Hixson-Crowell": "hixson_crowell",
    "Korsmeyer-Peppas_full": "korsmeyer_peppas_full",
}

# Golden parameter name -> engine parameter name
GOLDEN_PARAM_KEYS: dict[str, str] = {
    "k": "k1",
}

GOF_GOLDEN_KEYS = {"R2", "AIC", "RMSE", "MSC", "adj_R2"}

GOLDEN_QUANTITY_PREFIX = POLLI_SLOW_SIMILAR.label


def is_polli_golden(golden: dict[str, Any]) -> bool:
    return golden_section_key(golden) is not None


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
    status = "PASS" if abs_diff <= tol else "FAIL"
    return {
        "quantity": quantity,
        "computed": computed,
        "golden": golden_val,
        "abs_diff": abs_diff,
        "status": status,
    }


def _gof_from_summary(summary: ModelFitSummary, key: str) -> float | None:
    return summary.gof_mean.get(key)


def compare_polli_slow_similar_pair(
    pair: PairExportData,
    golden: dict[str, Any],
    *,
    default_tol: float = 0.01,
) -> list[dict[str, Any]]:
    """Compare polli_slow_similar against golden_polli_slow_similar.json only."""
    if not is_golden_polli_pair(pair.ref_file, pair.test_file):
        return []

    section_key = golden_section_key(golden)
    if section_key is None:
        return []

    tol = float(golden.get("_meta", {}).get("tolerance_deterministic", default_tol))
    section = golden[section_key]
    prefix = GOLDEN_QUANTITY_PREFIX
    rows: list[dict[str, Any]] = []

    for model_name, expected in section.get("reference_model_fits", {}).items():
        engine_key = GOLDEN_MODEL_KEYS.get(model_name)
        if engine_key is None:
            continue
        summary = pair.ref_summaries.get(engine_key)
        if summary is None:
            for key, golden_val in expected.items():
                rows.append(
                    _compare_value(
                        f"{prefix}.ref.{model_name}.{key}",
                        None,
                        float(golden_val),
                        tol,
                    )
                )
            continue

        for key, golden_val in expected.items():
            quantity = f"{prefix}.ref.{model_name}.{key}"
            if key in GOF_GOLDEN_KEYS:
                computed = _gof_from_summary(summary, key)
            else:
                engine_param = GOLDEN_PARAM_KEYS.get(key, key)
                computed = summary.param_mean.get(engine_param)
            rows.append(_compare_value(quantity, computed, float(golden_val), tol))

    sim = section.get("similarity_ref_vs_test_unit6", {})
    metrics = pair.metrics_fda
    for key in ("f2", "f1", "DE_ref", "DE_test", "MDT_ref", "MDT_test"):
        if key not in sim:
            continue
        rows.append(
            _compare_value(
                f"{prefix}.similarity.{key}",
                float(metrics[key]),
                float(sim[key]),
                tol,
            )
        )

    expected_times = sim.get("truncation_points_min")
    if expected_times is not None:
        computed_times = metrics.get("truncated_times", [])
        match = list(computed_times) == list(expected_times)
        rows.append(
            {
                "quantity": f"{prefix}.similarity.truncation_points_min",
                "computed": str(computed_times),
                "golden": str(expected_times),
                "abs_diff": None if match else float("nan"),
                "status": "PASS" if match else "FAIL",
            }
        )

    return rows


def compare_golden(
    pairs: list[PairExportData],
    golden: dict[str, Any],
    *,
    default_tol: float = 0.01,
) -> tuple[list[dict[str, Any]], bool]:
    """Run golden checks supported by the loaded golden file."""
    rows: list[dict[str, Any]] = []
    if is_polli_golden(golden):
        for pair in pairs:
            rows.extend(compare_polli_slow_similar_pair(pair, golden, default_tol=default_tol))
    else:
        flat = {k: v for k, v in golden.items() if not k.startswith("_")}
        tol = default_tol
        for quantity, golden_val in sorted(flat.items()):
            rows.append(_compare_value(quantity, None, float(golden_val), tol))

    all_pass = all(row["status"] == "PASS" for row in rows)
    return rows, all_pass
