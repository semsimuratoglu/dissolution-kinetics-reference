"""Extract golden expected values from a benchmark run."""

from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path


def main() -> None:
    root = Path(__file__).resolve().parent.parent
    env = {"PYTHONPATH": str(root)}
    cmd = [
        sys.executable,
        str(root / "run_benchmark.py"),
        "--data",
        str(root / "data"),
        "--pairs",
        "weibull_ref.csv:weibull_test.csv",
        "--bootstrap-iterations",
        "1000",
        "--seed",
        "42",
    ]
    proc = subprocess.run(cmd, capture_output=True, text=True, env={**subprocess.os.environ, **env})
    if proc.returncode != 0:
        print(proc.stderr or proc.stdout)
        raise SystemExit(proc.returncode)

    # Import after subprocess path is validated
    sys.path.insert(0, str(root))
    from refeng.io import load_dissolution_csv
    from refeng.metrics import compute_profile_metrics
    from refeng.models import fit_model_per_unit

    times, ref_units = load_dissolution_csv(root / "data" / "weibull_ref.csv")
    _, test_units = load_dissolution_csv(root / "data" / "weibull_test.csv")

    weibull_ref = fit_model_per_unit(times, ref_units, "weibull")
    metrics = compute_profile_metrics(times, ref_units, test_units, truncation="FDA")

    golden = {
        "weibull_ref_vs_weibull_test.ref.weibull.alpha_mean": weibull_ref.param_mean["alpha"],
        "weibull_ref_vs_weibull_test.ref.weibull.beta_mean": weibull_ref.param_mean["beta"],
        "weibull_ref_vs_weibull_test.ref.weibull.R2": weibull_ref.gof_mean["R2"],
        "weibull_ref_vs_weibull_test.ref.weibull.AIC": weibull_ref.gof_mean["AIC"],
        "weibull_ref_vs_weibull_test.f2": metrics["f2"],
        "weibull_ref_vs_weibull_test.f1": metrics["f1"],
        "weibull_ref_vs_weibull_test.DE_ref": metrics["DE_ref"],
        "weibull_ref_vs_weibull_test.DE_test": metrics["DE_test"],
        "weibull_ref_vs_weibull_test.MDT_ref": metrics["MDT_ref"],
        "weibull_ref_vs_weibull_test.MDT_test": metrics["MDT_test"],
    }

    out = root / "golden" / "expected.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w") as fh:
        json.dump(golden, fh, indent=2)
        fh.write("\n")

    print(json.dumps(golden, indent=2))
    print(f"\nWrote {out}")


if __name__ == "__main__":
    main()
