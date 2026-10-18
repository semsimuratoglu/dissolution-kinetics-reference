#!/usr/bin/env python3
"""Build the platform-vs-reference equivalence tables (CSV) from committed outputs only.

Inputs (no engine is re-run):
  platform_outputs/platform_results.json        platform export (7 pairs)
  results/similarity.csv, exposure.csv,
          model_fits.csv                        dual-engine reference tables (R + Python)
  results/reference_results.json                input sha256 per dataset; R values (bootstrap)
  results/python_reference_export.json          Python FDA-pathway f2/f1 and bootstrap

Outputs (CSV) in --out-dir:
  overview.csv, equivalence.csv, by_metric_class.csv,
  fda_pathway.csv, model_selection.csv, bootstrap_mendyk.csv

Run:
  PYTHONPATH=. python3 scripts/compare/build_equivalence_tables.py --out-dir /tmp/equivalence

Verdict per row from max(|Platform-R|, |Platform-Py|):
  <= 1e-6 epsilon; < 0.01 solver-level; >= 0.01 exceeds; missing Platform, R or Py value = n/a.
Values not present in an export are written as MISSING; nothing is filled by hand.
Not reproducible from public files (inputs not published): input-file identity against the
platform's own CSVs, and the platform MDT value from before the MDT fix.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import statistics
from pathlib import Path
from typing import Any, Optional

ROOT = Path(__file__).resolve().parents[2]

EPSILON = 1e-6
CRITERION = 0.01
MISSING = "MISSING"

# Dataset label -> platform export pair key.
PLATFORM_PAIR_KEYS = {
    "polli_slow_similar": "polli_slow_similar",
    "polli_fast_vs_slow_dissimilar": "polli_fast_vs_slow_dissimilar",
    "mendyk": "mendyk_reference_N12.csv:mendyk_test_N12.csv",
    "shah_similar": "shah_ref.csv:shah_test.csv",
    "shah_borderline": "shah_ref.csv:shah_borderline_test.csv",
    "tsong": "03a_tsong_REFERENCE.csv:03b_tsong_TEST.csv",
    "chow": "04a_chow_REFERENCE.csv:04b_chow_TEST.csv",
}

# (reference model name, platform model key, [(metric label, platform param key, reference param column)])
MODELS = [
    ("Zero-order", "Zero-order", [("k0", "k0", "k0")]),
    ("First-order", "First-order", [("k", "k", "k1")]),
    ("Higuchi", "Higuchi", [("kH", "kH", "kH")]),
    ("Hixson-Crowell", "Hixson-Crowell", [("kHC", "kHC", "kHC")]),
    ("Korsmeyer-Peppas (full)", "Korsmeyer-Peppas (full profile)", [("kKP", "kKP", "kKP"), ("n", "n", "n")]),
    ("Korsmeyer-Peppas (≤60%)", "Korsmeyer-Peppas", [("kKP", "kKP", "kKP"), ("n", "n", "n")]),
    ("Weibull", "Weibull", [("alpha", "alpha", "alpha"), ("beta", "beta", "beta")]),
]
GOF_METRICS = ["R2", "adjR2", "RMSE", "AIC"]

METRIC_CLASSES = [
    ("Similarity (f2, f1, CV)", "similarity"),
    ("Exposure (DE, MDT)", "exposure"),
    ("Goodness-of-fit (R2, adjR2, RMSE, AIC)", "gof"),
    ("Model parameters", "parameter"),
]

VERDICTS = ["epsilon", "solver-level", "exceeds", "n/a"]

# Python export keys (polli_fast_vs_slow_dissimilar is not exported).
PY_EXPORT_KEYS = {
    "polli_slow_similar": "polli_slow",
    "mendyk": "mendyk",
    "shah_similar": "shah",
    "shah_borderline": "shah_borderline",
    "tsong": "tsong",
    "chow": "chow",
}

PLATFORM_CANDIDATE_ALIASES = {"K-P full": "Korsmeyer-Peppas (full)"}
KP_LE60 = "Korsmeyer-Peppas (≤60%)"
DELTA_AIC = 2.0


def _num(value: Any) -> Optional[float]:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        return None
    return float(value)


def _csv_value(text: str) -> Any:
    # "{}" is the CSV writer's placeholder for an empty (non-converged) fit block.
    if text in ("", "{}"):
        return None
    try:
        return float(text)
    except ValueError:
        return text


def _read_table(path: Path) -> list[dict[str, Any]]:
    with path.open(newline="", encoding="utf-8") as fh:
        return [{k: _csv_value(v) for k, v in row.items()} for row in csv.DictReader(fh)]


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _diff(a: Optional[float], b: Optional[float]) -> Optional[float]:
    if a is None or b is None:
        return None
    return abs(a - b)


def _max(values: list[Optional[float]]) -> Optional[float]:
    vals = [v for v in values if v is not None]
    return max(vals) if vals else None


def _cell(value: Any) -> Any:
    return MISSING if value is None else value


def _write_csv(path: Path, rows: list[list[Any]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as fh:
        csv.writer(fh).writerows(rows)


def build_rows(platform: dict[str, Any], results: Path, ref_json: dict[str, Any]) -> list[dict[str, Any]]:
    similarity = {r["dataset"]: r for r in _read_table(results / "similarity.csv")}
    exposure = {r["dataset"]: r for r in _read_table(results / "exposure.csv")}
    fits = {(r["dataset"], r["model"]): r for r in _read_table(results / "model_fits.csv") if r["profile"] == "ref"}
    hashes = {d["label"]: (d["sha256_ref"], d["sha256_test"]) for d in ref_json["datasets"]}

    rows: list[dict[str, Any]] = []

    def add(dataset, arm, cls, metric, plat_v, r_v, py_v, input_key):
        rows.append(dict(dataset=dataset, arm=arm, cls=cls, metric=metric,
                         plat=_num(plat_v), r=_num(r_v), py=_num(py_v), input_key=input_key))

    for dataset, pair_key in PLATFORM_PAIR_KEYS.items():
        sp = platform["pairs"][pair_key]
        ref_hash, test_hash = hashes[dataset]
        pair_input = f"{ref_hash}|{test_hash}"
        sim = similarity[dataset]
        add(dataset, "pair", "similarity", "f2 (primary)", sp["similarity"]["f2_primary"],
            sim["f2_primary[R]"], sim["f2_primary[Py]"], pair_input)
        add(dataset, "pair", "similarity", "f1", sp["similarity"]["f1"],
            sim["f1_primary[R]"], sim["f1_primary[Py]"], pair_input)
        add(dataset, "pair", "similarity", "first_point_CV (%)", sp["cv_gate"]["first_point_cv"],
            sim["first_point_CV[R]"], sim["first_point_CV[Py]"], pair_input)

    for dataset, pair_key in PLATFORM_PAIR_KEYS.items():
        sp = platform["pairs"][pair_key]
        ref_hash, test_hash = hashes[dataset]
        exp = exposure[dataset]
        for metric in ("DE", "MDT"):
            for arm, arm_hash in (("ref", ref_hash), ("test", test_hash)):
                key = f"{metric}_{arm}"
                add(dataset, arm, "exposure", key, sp["exposure"][key],
                    exp[f"{key}[R]"], exp[f"{key}[Py]"], arm_hash)

    for dataset, pair_key in PLATFORM_PAIR_KEYS.items():
        sp = platform["pairs"][pair_key]
        ref_hash, _ = hashes[dataset]
        for ref_model, plat_model, params in MODELS:
            plat_fit = sp["reference_model_fits"].get(plat_model) or {}
            ref_fit = fits.get((dataset, ref_model), {})
            for gof in GOF_METRICS:
                add(dataset, ref_model, "gof", gof, plat_fit.get(gof),
                    ref_fit.get(f"{gof}[R]"), ref_fit.get(f"{gof}[Py]"), ref_hash)
            for label, plat_key, ref_col in params:
                add(dataset, ref_model, "parameter", label, plat_fit.get(plat_key),
                    ref_fit.get(f"{ref_col}[R]"), ref_fit.get(f"{ref_col}[Py]"), ref_hash)

    seen: set[tuple[str, str, str, str]] = set()
    for row in rows:
        row["d_plat_r"] = _diff(row["plat"], row["r"])
        row["d_plat_py"] = _diff(row["plat"], row["py"])
        row["d_r_py"] = _diff(row["r"], row["py"])
        if row["plat"] is None or row["r"] is None or row["py"] is None:
            row["max_plat_ref"] = None
            row["verdict"] = "n/a"
        else:
            m = max(row["d_plat_r"], row["d_plat_py"])
            row["max_plat_ref"] = m
            if m <= EPSILON:
                row["verdict"] = "epsilon"
            elif m < CRITERION:
                row["verdict"] = "solver-level"
            else:
                row["verdict"] = "exceeds"
        # Exposure metrics depend on one arm's file only, so DE_ref and DE_test share an identity.
        metric = row["metric"].split("_")[0] if row["cls"] == "exposure" else row["metric"]
        model = row["arm"] if row["cls"] in ("gof", "parameter") else ""
        identity = (row["cls"], model, metric, row["input_key"])
        row["shared_input_duplicate"] = identity in seen
        seen.add(identity)
    return rows


def write_core_tables(rows: list[dict[str, Any]], out_dir: Path, inputs: list[tuple[str, Path]]) -> dict[str, Any]:
    counts = {v: sum(r["verdict"] == v for r in rows) for v in VERDICTS}
    unique = {v: sum(r["verdict"] == v and not r["shared_input_duplicate"] for r in rows) for v in VERDICTS}
    dup_total = sum(r["shared_input_duplicate"] and r["verdict"] != "n/a" for r in rows)
    dup_all = sum(r["shared_input_duplicate"] for r in rows)
    gof = [r for r in rows if r["cls"] == "gof"]
    maxima = {
        "all": (_max([r["d_plat_r"] for r in rows]), _max([r["d_plat_py"] for r in rows]), _max([r["d_r_py"] for r in rows])),
        "gof": (_max([r["d_plat_r"] for r in gof]), _max([r["d_plat_py"] for r in gof]), _max([r["d_r_py"] for r in gof])),
    }

    ov = [["Platform vs independent reference (R + Python) — equivalence table"],
          ["Values read from committed outputs; no engine was re-run."], [],
          ["Input", "Path", "sha256"]]
    ov += [[label, str(path.relative_to(ROOT)) if path.is_relative_to(ROOT) else path.name, _sha256(path)]
           for label, path in inputs]
    ov += [[], ["Verdict rule", "max(|Platform-R|, |Platform-Py|): <= 1e-6 epsilon; < 0.01 solver-level; >= 0.01 exceeds; missing Platform, R or Py = n/a"],
           [], ["Verdict", "All rows", "Unique rows (shared_input_duplicate = FALSE)"]]
    ov += [[v, counts[v], unique[v]] for v in VERDICTS]
    ov += [["Total", len(rows), sum(unique.values())],
           ["shared_input_duplicate (non-n/a rows)", dup_total],
           ["shared_input_duplicate (all rows)", dup_all], [],
           ["Maxima", "max |Platform-R|", "max |Platform-Py|", "max |R-Py|"],
           ["All quantities", *maxima["all"]],
           ["Goodness-of-fit only", *maxima["gof"]]]
    _write_csv(out_dir / "overview.csv", ov)

    eq = [["Dataset", "Arm / Model", "Class", "Metric", "Platform", "R", "Py", "|Platform-R|", "|Platform-Py|",
           "|R-Py|", "max(|Platform-R|,|Platform-Py|)", "Verdict", "shared_input_duplicate", "Input sha256"]]
    eq += [[r["dataset"], r["arm"], r["cls"], r["metric"], r["plat"], r["r"], r["py"], r["d_plat_r"],
            r["d_plat_py"], r["d_r_py"], r["max_plat_ref"], r["verdict"], r["shared_input_duplicate"],
            r["input_key"]] for r in rows]
    _write_csv(out_dir / "equivalence.csv", eq)

    cls_rows = [["Metric class", "n (non-n/a)", "median max(|Platform-R|,|Platform-Py|)",
                 "max max(|Platform-R|,|Platform-Py|)", "max |R-Py|", "count >= 0.01"]]
    class_summary = []
    for label, key in METRIC_CLASSES:
        members = [r for r in rows if r["cls"] == key and r["verdict"] != "n/a"]
        dev = [r["max_plat_ref"] for r in members]
        entry = [label, len(members), statistics.median(dev) if dev else None, max(dev) if dev else None,
                 _max([r["d_r_py"] for r in members]), sum(d >= CRITERION for d in dev)]
        class_summary.append(entry)
        cls_rows.append(entry)
    _write_csv(out_dir / "by_metric_class.csv", cls_rows)
    return {"counts": counts, "unique": unique, "duplicates_non_na": dup_total,
            "duplicates_all": dup_all, "maxima": maxima, "classes": class_summary}


def write_fda_pathway(platform: dict[str, Any], py_export: dict[str, Any], ref_json: dict[str, Any], out_dir: Path) -> list[list[Any]]:
    r_block = ref_json["r"]
    rows = [["Dataset", "Metric", "Platform", "R", "Py", "|Platform-R|", "|Platform-Py|", "|R-Py|",
             "Platform FDA retained times", "Note"]]
    for dataset, pair_key in PLATFORM_PAIR_KEYS.items():
        sim = platform["pairs"][pair_key]["similarity"]
        py_key = PY_EXPORT_KEYS.get(dataset)
        py_fda = py_export.get(py_key, {}).get("fda", {}) if py_key else {}
        r_sim = r_block.get(dataset, {}).get("similarity", {})
        for metric in ("f2", "f1"):
            p, r, y = _num(sim.get(f"{metric}_fda")), _num(r_sim.get(f"{metric}_fda")), _num(py_fda.get(metric))
            notes = []
            if p is None:
                notes.append("platform export has no FDA-pathway " + metric)
            if r is None:
                notes.append("R export has no FDA-pathway fields")
            if y is None:
                notes.append("pair not in Python export" if py_key is None else "Python export has no FDA " + metric)
            rows.append([dataset, f"{metric} (FDA)", _cell(p), _cell(r), _cell(y), _cell(_diff(p, r)),
                         _cell(_diff(p, y)), _cell(_diff(r, y)),
                         ", ".join(f"{t:g}" for t in sim.get("trunc_fda_min", [])), "; ".join(notes)])
    _write_csv(out_dir / "fda_pathway.csv", rows)
    return rows[1:]


def _candidate_set(aics: dict[str, float]) -> list[str]:
    if not aics:
        return []
    best = min(aics.values())
    return sorted((m for m, a in aics.items() if a - best < DELTA_AIC), key=lambda m: aics[m])


def write_model_selection(platform: dict[str, Any], results: Path, out_dir: Path) -> list[list[Any]]:
    fits = {(r["dataset"], r["model"]): r for r in _read_table(results / "model_fits.csv") if r["profile"] == "ref"}
    plat_to_ref = {plat_model: ref_model for ref_model, plat_model, _ in MODELS}
    rows = [["Dataset", "Engine", "N", "Models considered", "Candidate set (by AIC)", "AIC of candidates",
             "Platform exported candidate set", "Recomputed == platform exported", "Identical across engines"]]
    for dataset, pair_key in PLATFORM_PAIR_KEYS.items():
        sp = platform["pairs"][pair_key]
        ranking = sp.get("reference_aic_ranking") or {}
        fit_n = ranking.get("reference_fit_n")
        plat_aics: dict[str, float] = {}
        for ref_model, plat_model, _ in MODELS:
            if ref_model == KP_LE60:
                continue
            f = sp["reference_model_fits"].get(plat_model) or {}
            if _num(f.get("AIC")) is not None and f.get("fitN") == fit_n:
                plat_aics[ref_model] = float(f["AIC"])
        engines = {"Platform": plat_aics}
        for eng, col in (("R", "AIC[R]"), ("Py", "AIC[Py]")):
            aics = {}
            for ref_model, _, _ in MODELS:
                if ref_model == KP_LE60:
                    continue
                v = _num(fits.get((dataset, ref_model), {}).get(col))
                if v is not None:
                    aics[ref_model] = v
            engines[eng] = aics
        sets = {eng: _candidate_set(a) for eng, a in engines.items()}
        identical = len({tuple(sorted(s)) for s in sets.values()}) == 1
        exported = [PLATFORM_CANDIDATE_ALIASES.get(c["model"], plat_to_ref.get(c["model"], c["model"]))
                    for c in ranking.get("candidates", [])]
        for eng in ("Platform", "R", "Py"):
            cand = sets[eng]
            rows.append([dataset, eng, fit_n if eng == "Platform" else "not exported", len(engines[eng]),
                         ", ".join(cand), ", ".join(f"{engines[eng][m]:.6f}" for m in cand),
                         ", ".join(exported) if eng == "Platform" else "",
                         (sorted(cand) == sorted(exported)) if eng == "Platform" else "", identical])
    _write_csv(out_dir / "model_selection.csv", rows)
    return rows[1:]


def write_bootstrap_mendyk(platform: dict[str, Any], py_export: dict[str, Any], ref_json: dict[str, Any], out_dir: Path) -> list[list[Any]]:
    pb_ = platform["pairs"][PLATFORM_PAIR_KEYS["mendyk"]].get("bootstrap") or {}
    py_b = py_export.get("mendyk", {}).get("boot_ema", {})
    r_b = ref_json["r"].get("mendyk", {}).get("bootstrap") or {}
    rows = [["Quantity", "Platform", "R", "Py", "|Platform-R|", "|Platform-Py|", "|R-Py|"]]
    for label, key in (("E_f2", "E_f2"), ("CI lower", "ci_lower"), ("CI upper", "ci_upper"),
                       ("bc_f2", "bc_f2"), ("median", "median")):
        p, r, y = _num(pb_.get(key)), _num(r_b.get(key)), _num(py_b.get(key))
        rows.append([label, _cell(p), _cell(r), _cell(y), _cell(_diff(p, r)), _cell(_diff(p, y)), _cell(_diff(r, y))])
    _write_csv(out_dir / "bootstrap_mendyk.csv", rows)
    return rows[1:]


def main(argv: Optional[list[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--platform", default=str(ROOT / "platform_outputs" / "platform_results.json"))
    parser.add_argument("--results", default=str(ROOT / "results"), help="Folder with the reference result tables")
    parser.add_argument("--out-dir", required=True, help="Output folder for the CSV tables")
    args = parser.parse_args(argv)

    platform_path = Path(args.platform).resolve()
    results = Path(args.results).resolve()
    out_dir = Path(args.out_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    ref_json_path = results / "reference_results.json"
    py_path = results / "python_reference_export.json"

    platform = json.loads(platform_path.read_text())
    ref_json = json.loads(ref_json_path.read_text())
    py_export = json.loads(py_path.read_text())

    rows = build_rows(platform, results, ref_json)
    inputs = [("platform", platform_path)] + [(p.name, p) for p in (
        results / "similarity.csv", results / "exposure.csv", results / "model_fits.csv", ref_json_path, py_path)]
    summary = write_core_tables(rows, out_dir, inputs)

    print(f"Rows: {len(rows)}")
    for v in VERDICTS:
        print(f"  {v:<13} all={summary['counts'][v]:<4} unique={summary['unique'][v]}")
    print(f"  shared_input_duplicate: {summary['duplicates_non_na']} (non-n/a), {summary['duplicates_all']} (all rows)")
    for scope, (a, b, c) in summary["maxima"].items():
        print(f"  max [{scope}]: |Platform-R|={a:.4g} |Platform-Py|={b:.4g} |R-Py|={c:.4g}")
    for label, n, med, mx, rpy, exceed in summary["classes"]:
        print(f"  {label}: n={n} median={med:.3g} max={mx:.4g} max|R-Py|={rpy:.4g} >=0.01: {exceed}")

    print("FDA_Pathway:")
    for r in write_fda_pathway(platform, py_export, ref_json, out_dir):
        print(f"  {r[0]:<30} {r[1]:<9} Platform={r[2]} R={r[3]} Py={r[4]} |Platform-Py|={r[6]}")
    print("Model_Selection:")
    for r in write_model_selection(platform, results, out_dir):
        extra = f" exported=[{r[6]}] match={r[7]}" if r[1] == "Platform" else ""
        print(f"  {r[0]:<30} {r[1]:<8} [{r[4]}] identical_across_engines={r[8]}{extra}")
    print("Bootstrap_mendyk:")
    for r in write_bootstrap_mendyk(platform, py_export, ref_json, out_dir):
        print(f"  {r[0]:<9} Platform={r[1]} R={r[2]} Py={r[3]} |Platform-R|={r[4]} |Platform-Py|={r[5]} |R-Py|={r[6]}")
    print(f"Wrote CSV tables to: {out_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
