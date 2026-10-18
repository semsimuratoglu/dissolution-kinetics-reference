# Independent Reference Engines for Dissolution Kinetics and f2 Similarity (Python & R)

Independent reference implementations for dissolution profile kinetic modelling,
similarity metrics, and bootstrap uncertainty quantification. The Python and R
engines were implemented from first principles using standard scientific libraries
only.

**Licenses:** source code is [BSD-3-Clause](LICENSE); datasets under `data/` and
derived result tables under `results/` are [CC-BY-4.0](https://creativecommons.org/licenses/by/4.0/).

## What the engines compute

### Kinetic model fitting

Each experimental unit (replicate row in the CSV) is fitted independently. Reported
parameters are the mean ± SD across units. Fits use residuals on **Q (percent dissolved)**,
never a linearised transform.

| Model | Equation (Q in %) |
|-------|-------------------|
| Zero-order | Q = k₀ · t |
| First-order | Q = 100 · (1 − e^(−k₁ · t)) |
| Higuchi | Q = k_H · √t |
| Hixson–Crowell | 100^(1/3) − (100 − Q)^(1/3) = k_HC · t |
| Korsmeyer–Peppas | Q = k_KP · t^n (full profile or ≤60% window) |
| Weibull | Q = 100 · (1 − exp(−(t/α)^β)) |

Unbounded fits use Levenberg–Marquardt (`method='lm'`); positivity-constrained
parameters use `method='trf'` with bounds.

### Goodness-of-fit

Per model and arm: R², adjusted R², RMSE = √(RSS/(n−p)), AIC = N·ln(RSS)+2P
(DDSolver convention), and MSC.

### Similarity and exposure metrics

- **f₂** and **f₁** with FDA and EMA 85% truncation rules:
  - FDA: retain points up to the first time where **both** reference and test mean exceed 85%.
  - EMA: retain points up to and including the first time where **either** reference or test mean exceeds 85%
    (EMA BE guideline Appendix I: not more than one mean value > 85% for either formulation).
- **DE%** (dissolution efficiency) and **MDT** (mean dissolution time) on the full profile.

### Bootstrap f₂

When the %CV gate is exceeded (first time-point CV > 20% or any later point > 10%),
unit-level bootstrap resampling (10 000 iterations, seed 42) yields E(f₂), bc-f₂, and a
90% percentile confidence interval. Similarity when the lower 90% CI bound ≥ 50.

## Requirements

- **Python 3.9+** (tested with 3.9.6) (core engine)
- **R 4.x** with `minpack.lm` and `bootf2` (optional independent cross-check)

## Quick start (Python)

```bash
python -m pip install -r requirements.txt
PYTHONPATH=. python run_benchmark.py \
  --data data \
  --pairs "polli_slow_similar_ref.csv:polli_slow_similar_test.csv,mendyk_reference_N12.csv:mendyk_test_N12.csv,shah_ref.csv:shah_test.csv" \
  --golden golden/golden_polli_slow_similar.json \
  --golden-shah-external golden/golden_shah_external.json
```

Run the test suite:

```bash
PYTHONPATH=. python -m pytest tests/ -v
```

Regenerate dual-engine result exports (JSON + CSV):

```bash
PYTHONPATH=. python scripts/export_reference_results.py
PYTHONPATH=. python scripts/export_python_reference.py
PYTHONPATH=. python scripts/generate_reference_export.py
```

## Quick start (R, optional)

R is not required for the Python workflow. When Rscript is available:

```bash
bash scripts/run_r_crosscheck.sh data "polli_slow_similar_ref.csv:polli_slow_similar_test.csv"
# or
PYTHONPATH=. python run_benchmark.py ... --with-r
```

Install R packages: `install.packages(c('minpack.lm','bootf2','jsonlite'))`.

Export Shah 1998 datasets from bootf2: `Rscript scripts/export_shah1998.R --out-dir data`.

## Project layout

| Path | Purpose |
|------|---------|
| `refeng/` | Python reference engine modules |
| `run_benchmark.py` | CLI and golden-value regression |
| `reference.R` | Independent R engine (`nlsLM` + `bootf2`) |
| `data/` | Input CSV fixtures (`Time`, `Rep 1`…`Rep N`) |
| `golden/` | Externally anchored expected-value JSON |
| `results/` | Committed dual-engine JSON/CSV exports |
| `scripts/` | Export and data-generation utilities |
| `dev/` | Development-only validation reports (not required for reproduction) |

See `data/DATA_README.md` for dataset provenance, including the three distinct
Polli comparison sets.

## Dataset provenance

| Set | Source | Notes |
|-----|--------|-------|
| polli_slow_similar | Polli 1997; Zhang 2010 DDSolver Table 2 | Golden validation pair (f₂ ≈ 91) |
| polli_fast_hc_demo | Polli benchmark | Reference-only HC demo |
| polli_fast_vs_slow_dissimilar | Polli 1997 benchmark | Dissimilarity demo (f₂ ≈ 22) |
| shah_* | bootf2::shah1998 v0.4.1 | Literature-anchored f₂ values |
| tsong | Tsong 1996 | R=Python cross-check |
| chow | Chow & Ki 1997 Table 1 | Hours-based profile; R=Python cross-check |
| **mendyk N12** | **Synthetic fixture** | **NOT** Mendyk 2013 published raw data; engineered high-CV pair (times 15–120 min) to exercise bootstrap |

## Pinned dependencies

**Python** (`requirements.txt`): numpy 2.0.2, scipy 1.13.1, pandas 2.2.3, pytest 8.3.4, openpyxl 3.1.5.

**R** (`renv-versions.txt`): minpack.lm 1.2-4, bootf2 0.4.1, jsonlite 2.0.0.

## Citation

See [CITATION.cff](CITATION.cff) for metadata suitable for Zenodo and reference managers.
