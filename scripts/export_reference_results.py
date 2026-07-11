#!/usr/bin/env python3
"""Generate dual-engine reference results as JSON and CSV exports."""

from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from refeng.dual_export import write_reference_results_exports  # noqa: E402


def main() -> int:
    data_dir = ROOT / "data"
    output_dir = ROOT / "results"
    script_dir = ROOT / "scripts"

    paths = write_reference_results_exports(output_dir, data_dir, script_dir)
    for name, path in paths.items():
        print(f"{name}: {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
