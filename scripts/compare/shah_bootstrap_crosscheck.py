#!/usr/bin/env python3
"""Shah 1998 bootstrap cross-check: platform output vs Python refeng, R and bootf2 (CSV output).

Inputs:
  platform_outputs/platform_shah_bootstrap.json   platform bootstrap export (seed 42 + seeds 42-66)
  data/                                           Shah CSVs (byte-identity with the platform inputs is checked)

Engines (all seed 42, 10,000 iterations, 90% percentile CI):
  Platform  read from the platform export
  Python    refeng.bootstrap.bootstrap_f2(truncation="EMA")
  R         reference.R python_style_bootstrap (only function definitions are evaluated)
  bootf2    bootf2::bootf2(f2.type="exp.f2", ci.type="percentile", quantile.type="7"):
              regulation="EMA"                       -> decision context
              regulation="FDA", both.TR.85=TRUE       -> context only

Outputs (CSV) in --out-dir: readme.csv, engines.csv, seed_stability.csv, bootf2_configurations.csv

Run:
  PYTHONPATH=. python3 scripts/compare/shah_bootstrap_crosscheck.py --out-dir /tmp/shah_crosscheck
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import platform as py_platform
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import Any, Optional

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

import numpy as np  # noqa: E402

from refeng.bootstrap import bootstrap_f2  # noqa: E402
from refeng.io import load_dissolution_csv  # noqa: E402
from refeng.metrics import ema_truncation, fda_truncation  # noqa: E402

N_ITER = 10_000
SEED = 42

DATASETS = {
    "shah_similar": ("shah_ref.csv", "shah_test.csv"),
    "shah_borderline": ("shah_ref.csv", "shah_borderline_test.csv"),
}

# Predicted platform values from an independent port of the platform bootstrap function
# (same PRNG stream, EMA either-product cut). A prediction, not an oracle.
PREDICTION = {
    "shah_similar": {"E_f2": 56.8316, "ci_lower": 50.2413, "ci_upper": 66.6034},
    "shah_borderline": {"E_f2": 45.1635, "ci_lower": 43.2790, "ci_upper": 47.9423},
}

R_SCRIPT = r"""
suppressPackageStartupMessages(library(jsonlite))
a <- commandArgs(trailingOnly = TRUE)
reference_r <- a[1]; data_dir <- a[2]; n_iter <- as.integer(a[3]); seed <- as.integer(a[4])
labels <- strsplit(a[5], ",")[[1]]; refs <- strsplit(a[6], ",")[[1]]; tests <- strsplit(a[7], ",")[[1]]

for (e in parse(reference_r)) {
  if (is.call(e) && identical(e[[1]], as.name("<-")) && is.name(e[[2]])) {
    rhs <- e[[3]]
    if ((is.call(rhs) && identical(rhs[[1]], as.name("function"))) ||
        as.character(e[[2]]) %in% c("BOOT_ITER", "BOOT_SEED")) eval(e, globalenv())
  }
}

read_units <- function(f) {
  d <- read.csv(file.path(data_dir, f), check.names = FALSE)
  list(times = d[[1]], units = t(as.matrix(d[, -1, drop = FALSE])), df = d)
}

bootf2_run <- function(r, t, regulation, both) {
  rd <- r$df; td <- t$df
  names(rd)[1] <- "time"; names(td)[1] <- "time"
  res <- bootf2::bootf2(td, rd, n.boots = n_iter, seed = seed, regulation = regulation,
                        both.TR.85 = both, f2.type = "exp.f2", ci.type = "percentile",
                        quantile.type = "7", digits = 10, time.unit = "min",
                        print.report = FALSE, output.to.screen = FALSE)
  rm <- as.matrix(rd); tm <- as.matrix(td)
  cf <- bootf2:::calcf2(test = tm, ref = rm, regulation = regulation, digits = 10, cv.rule = FALSE,
                        min.points = 1, both.TR.85 = both, message = FALSE, f2.type = "exp.f2",
                        plot = FALSE, time.unit = "min")
  s <- res$boot.summary[res$boot.summary$f2.type == "exp.f2", ]
  ci <- res$boot.ci[res$boot.ci$f2.type == "exp.f2" & grepl("Percentile", res$boot.ci$ci.type), ]
  tp <- table(res$boot.f2$f2.tp)
  list(regulation = regulation, both_TR_85 = both, ci_type = ci$ci.type,
       E_f2 = s$f2o, boot_mean = s$boot.mean, boot_median = s$boot.median, boot_na = s$boot.na,
       ci_lower = ci$ci.lower, ci_upper = ci$ci.upper,
       points_kept = unname(cf$f2.tp), retained_times = r$times[seq_len(cf$f2.tp)],
       draw_points_distribution = as.list(setNames(as.integer(tp), names(tp))))
}

out <- list(r_version = R.version.string,
            bootf2_version = as.character(packageVersion("bootf2")),
            datasets = list())
for (i in seq_along(labels)) {
  r <- read_units(refs[i]); t <- read_units(tests[i])
  ema <- ema_truncate(r$times, colMeans(r$units), colMeans(t$units))
  fda <- fda_truncate(r$times, colMeans(r$units), colMeans(t$units))
  pb <- python_style_bootstrap(r$times, r$units, t$units, n_iter = n_iter, seed = seed)
  out$datasets[[labels[i]]] <- list(
    python_style = c(pb, list(retained_times = ema$times,
                              primary_rule = if (identical(fda$idx, ema$idx)) "FDA" else "EMA")),
    bootf2_ema = bootf2_run(r, t, "EMA", FALSE),
    bootf2_fda = bootf2_run(r, t, "FDA", TRUE))
}
cat(toJSON(out, auto_unbox = TRUE, digits = NA, pretty = TRUE))
"""


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _git_head() -> Optional[str]:
    try:
        return subprocess.run(["git", "rev-parse", "HEAD"], cwd=ROOT, capture_output=True,
                              text=True, check=True).stdout.strip()
    except (OSError, subprocess.CalledProcessError):
        return None


def _decision(ci_lower: Optional[float]) -> str:
    if ci_lower is None:
        return "MISSING"
    return "similar" if ci_lower >= 50.0 else "not similar"


def _times(values) -> str:
    return ", ".join(f"{float(v):g}" for v in values)


def _diff(a: Optional[float], b: Optional[float]) -> Optional[float]:
    return None if a is None or b is None else abs(a - b)


def _write_csv(path: Path, rows: list[list[Any]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as fh:
        csv.writer(fh).writerows(rows)


def run_python(data: Path) -> dict[str, Any]:
    out = {}
    for label, (ref_file, test_file) in DATASETS.items():
        times, ref_units = load_dissolution_csv(data / ref_file)
        _, test_units = load_dissolution_csv(data / test_file)
        ref_means, test_means = ref_units.mean(axis=0), test_units.mean(axis=0)
        ema = ema_truncation(times, ref_means, test_means)
        fda = fda_truncation(times, ref_means, test_means)
        boot = bootstrap_f2(times, ref_units, test_units, n_iterations=N_ITER, seed=SEED, truncation="EMA")
        out[label] = {
            "retained_times": ema.times.tolist(),
            "primary_rule": "FDA" if np.array_equal(fda.kept_indices, ema.kept_indices) else "EMA",
            "f2_observed": boot.f2_observed, "E_f2": boot.e_f2, "boot_mean": boot.mean,
            "median": boot.median, "ci_lower": boot.ci_lower, "ci_upper": boot.ci_upper,
            "bc_f2": boot.bc_f2, "truncation": boot.truncation,
        }
    return out


def run_r(data: Path) -> dict[str, Any]:
    with tempfile.TemporaryDirectory() as tmp:
        script = Path(tmp) / "shah_crosscheck.R"
        script.write_text(R_SCRIPT)
        labels = list(DATASETS)
        proc = subprocess.run(
            ["Rscript", str(script), str(ROOT / "reference.R"), str(data), str(N_ITER), str(SEED),
             ",".join(labels), ",".join(DATASETS[k][0] for k in labels),
             ",".join(DATASETS[k][1] for k in labels)],
            cwd=tmp, capture_output=True, text=True,
        )
    if proc.returncode != 0:
        raise RuntimeError(f"Rscript failed:\n{proc.stderr}")
    return json.loads(proc.stdout)


def write_tables(plat: dict[str, Any], plat_path: Path, py: dict[str, Any], r: dict[str, Any],
                 data: Path, out_dir: Path) -> list[list[Any]]:
    plat_by_label = {d["label"]: d for d in plat["datasets"]}
    seeds = plat["datasets"][0]["seed_sweep"]["seeds"]
    readme = [
        ["Shah 1998 bootstrap cross-check (shah_similar = ref vs test1, shah_borderline = ref vs test5)"],
        [],
        ["Iterations", N_ITER],
        ["Seed (Python, R, bootf2, platform default)", SEED],
        ["Platform seed sweep", f"{seeds[0]}–{seeds[-1]}"],
        ["CI", "90% percentile; decision: CI lower >= 50"],
        ["Platform", f"{plat['resampling']}; {plat['truncation']}"],
        ["Platform commit", plat.get("git_head")],
        ["Platform export", plat_path.name],
        ["Platform export sha256", _sha256(plat_path)],
        ["Platform export generated_at", plat.get("generated_at")],
        ["Reference git HEAD", _git_head()],
        ["Python", f"{py_platform.python_version()} (numpy {np.__version__})"],
        ["Python engine", "refeng.bootstrap.bootstrap_f2(truncation='EMA'); numpy default_rng"],
        ["R", r["r_version"]],
        ["R engine", "reference.R python_style_bootstrap; set.seed + sample()"],
        ["bootf2", r["bootf2_version"]],
        ["bootf2 call", "bootf2(test, ref, n.boots=10000, seed=42, f2.type='exp.f2', ci.type='percentile', quantile.type='7', digits=10)"],
        ["Prediction", "Platform E-f2 / CI predicted by an independent port of the platform function — a prediction, not an oracle"],
        [],
        ["Dataset", "Arm", "File", "sha256 (data/)", "sha256 (platform input)", "Identical"],
    ]
    for label, (ref_file, test_file) in DATASETS.items():
        pd_ = plat_by_label[label]
        for arm, f, plat_sha in (("ref", ref_file, pd_["ref_sha256"]), ("test", test_file, pd_["test_sha256"])):
            local = _sha256(data / f)
            readme.append([label, arm, f, local, plat_sha, local == plat_sha])
    _write_csv(out_dir / "readme.csv", readme)

    engines = [["Dataset", "Engine", "Retained times (observed)", "f2 observed", "E-f2", "CI lower", "CI upper",
                "Median", "bc-f2", "Decision", "Predicted E-f2", "Predicted CI lower", "Predicted CI upper",
                "|E-f2 - predicted|", "|CI lower - predicted|", "|CI upper - predicted|"]]
    for label in DATASETS:
        pd_ = plat_by_label[label]
        b = pd_["bootstrap_seed_default"]
        x = PREDICTION[label]
        ps = r["datasets"][label]["python_style"]
        be = r["datasets"][label]["bootf2_ema"]
        rows = [
            ["Platform", pd_["ema_retained_times"], b["f2_observed"], b["E_f2"], b["ci_lower"], b["ci_upper"],
             b["median"], b["bc_f2"]],
            ["Python refeng", py[label]["retained_times"], py[label]["f2_observed"], py[label]["E_f2"],
             py[label]["ci_lower"], py[label]["ci_upper"], py[label]["median"], py[label]["bc_f2"]],
            ["R python_style_bootstrap", ps["retained_times"], ps["f2_observed"], ps["E_f2"], ps["ci_lower"],
             ps["ci_upper"], ps["median"], ps["bc_f2"]],
            ["bootf2 EMA exp.f2", be["retained_times"], None, be["E_f2"], be["ci_lower"], be["ci_upper"],
             be["boot_median"], None],
        ]
        for name, times, f2o, ef2, lo, hi, med, bc in rows:
            is_plat = name == "Platform"
            engines.append([label, name, _times(times), f2o if f2o is not None else "n/a", ef2, lo, hi, med,
                            bc if bc is not None else "n/a", _decision(lo),
                            x["E_f2"] if is_plat else "", x["ci_lower"] if is_plat else "",
                            x["ci_upper"] if is_plat else "",
                            _diff(ef2, x["E_f2"]) if is_plat else "", _diff(lo, x["ci_lower"]) if is_plat else "",
                            _diff(hi, x["ci_upper"]) if is_plat else ""])
    _write_csv(out_dir / "engines.csv", engines)

    seed_rows = [["Dataset", "CI lower min", "CI lower max", "CI lower mean", "CI lower SD",
                  "Seeds with CI lower < 50", "Seeds"]]
    for label in DATASETS:
        s = plat_by_label[label]["seed_sweep"]
        cs = s["ci_lower_summary"]
        seed_rows.append([label, cs["min"], cs["max"], cs["mean"], cs["sd"], cs["n_below_50"], len(s["seeds"])])
    seed_rows += [[], ["Dataset", "Seed", "E-f2", "CI lower", "CI upper", "Median", "Decision"]]
    for label in DATASETS:
        for run in plat_by_label[label]["seed_sweep"]["runs"]:
            seed_rows.append([label, run["seed"], run["E_f2"], run["ci_lower"], run["ci_upper"], run["median"],
                              _decision(run["ci_lower"])])
    _write_csv(out_dir / "seed_stability.csv", seed_rows)

    cfg = [["Dataset", "Configuration", "Points kept (observed)", "Retained times (observed)", "E-f2 (f2o)",
            "Boot mean", "Boot median", "CI lower", "CI upper", "Decision", "Points used per draw (count)", "boot.na"]]
    for label in DATASETS:
        for key, name in (("bootf2_ema", "regulation=EMA (decision)"),
                          ("bootf2_fda", "regulation=FDA, both.TR.85=TRUE (context only)")):
            c = r["datasets"][label][key]
            dist = ", ".join(f"{k} pts: {v}" for k, v in c["draw_points_distribution"].items())
            cfg.append([label, name, c["points_kept"], _times(c["retained_times"]), c["E_f2"], c["boot_mean"],
                        c["boot_median"], c["ci_lower"], c["ci_upper"], _decision(c["ci_lower"]), dist, c["boot_na"]])
    _write_csv(out_dir / "bootf2_configurations.csv", cfg)
    return engines[1:]


def main(argv: Optional[list[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--platform", default=str(ROOT / "platform_outputs" / "platform_shah_bootstrap.json"))
    parser.add_argument("--data", default=str(ROOT / "data"), help="Reference CSV folder")
    parser.add_argument("--out-dir", required=True, help="Output folder for the CSV tables")
    args = parser.parse_args(argv)

    plat_path = Path(args.platform).resolve()
    data = Path(args.data).resolve()
    out_dir = Path(args.out_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    plat = json.loads(plat_path.read_text())

    py = run_python(data)
    r = run_r(data)
    rows = write_tables(plat, plat_path, py, r, data, out_dir)

    for row in rows:
        print(f"{row[0]:<16} {row[1]:<26} times=[{row[2]}] E-f2={row[4]:.4f} "
              f"CI={row[5]:.4f}–{row[6]:.4f} {row[9]}")
    for label in DATASETS:
        for key in ("bootf2_ema", "bootf2_fda"):
            c = r["datasets"][label][key]
            print(f"{label:<16} {key:<10} points={c['points_kept']} E-f2={c['E_f2']:.4f} "
                  f"CI={c['ci_lower']:.4f}–{c['ci_upper']:.4f} {_decision(c['ci_lower'])} "
                  f"draws={c['draw_points_distribution']}")
    print(f"Wrote CSV tables to: {out_dir}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
