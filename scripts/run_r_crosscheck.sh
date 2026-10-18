#!/usr/bin/env bash
# Optional R cross-check for the benchmark reference engine.
# Non-blocking: exits 0 when Rscript is missing (prints install hint).

set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
DATA_DIR="${1:-$ROOT/data}"
PAIRS="${2:-polli_slow_similar_ref.csv:polli_slow_similar_test.csv}"

if ! command -v Rscript >/dev/null 2>&1; then
  echo "R cross-check skipped (Rscript not found)."
  echo "Install R, then run: install.packages(c('minpack.lm','bootf2'))"
  exit 0
fi

Rscript "$ROOT/reference.R" --data "$DATA_DIR" --pairs "$PAIRS"
