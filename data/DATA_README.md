# Benchmark data for the independent reference engine

CSV format: `Time, Rep 1 … Rep N` (first column time; remaining columns are individual replicate units).

## Polli sets — three distinct comparisons (do not conflate)

| Label | Ref CSV | Test CSV | Publication | Profile | Role | f₂ (primary) |
|-------|---------|----------|-------------|---------|------|--------------|
| **polli_slow_similar** | `polli_slow_similar_ref.csv` | `polli_slow_similar_test.csv` (Unit 6 only) | Polli 1997; Zhang 2010 DDSolver Table 2 | Slow metoprolol — 6-unit mean vs slow Unit 6 | Golden / deterministic-equivalence validation | **91.02** SIMILAR |
| **polli_fast_hc_demo** | `polli_fast_hc_demo_ref.csv` | — (reference only) | Polli benchmark / Lopressor-like fast | Fast formulation (near-complete by 15 min) | Hixson–Crowell unconstrained demo (kHC≈0.033, R²≈0.83) | — |
| **polli_fast_vs_slow_dissimilar** | `polli_fast_vs_slow_dissimilar_ref.csv` | `polli_fast_vs_slow_dissimilar_test.csv` (6 slow units) | Polli 1997 benchmark | Fast ref vs slow test | Dissimilarity demonstration pair | **22.49** DISSIMILAR |

Legacy filenames (`polli_ref.csv`, `01_polli_*`, `02a_*`, `02b_*`) are **removed**. Code resolves
canonical labels via `refeng/polli_sets.py`.

## Other benchmark files

| File | Units | Time unit | Notes |
|---|---|---|---|
| tsong_ref.csv / tsong_test.csv | 6 | min | Tsong 1996. |
| chow_ref.csv / chow_test.csv | 12 | HOURS | Chow & Ki 1997 Table 1 (Tsong–Hammerstrom); times 1,2,3,4,6,8,10 h; 12 units each arm. |
| mendyk_reference_N12.csv / mendyk_test_N12.csv | 12 | min | **Synthetic high-variability set** (NOT Mendyk 2013 published raw data; times 15–120 min vs published 3–30 min). First-point CV ≈ 25.7% → forces bootstrap. |
| shah_ref.csv / shah_test.csv | 12 | min | **bootf2::shah1998** — `ref` vs `test1` (shah_similar); f₂≈60.03; bootstrap skipped. |
| shah_ref.csv / shah_borderline_test.csv | 12 | min | **bootf2::shah1998** — `ref` vs `test5` (shah_borderline); f₂≈48.05; bootstrap skipped. |

Export all Shah arms: `Rscript scripts/export_shah1998.R --out-dir data`.

## External Shah golden (`golden/golden_shah_external.json`)

Literature oracle from **Noce 2020 Table 5** and **Boddu 2024 Table IV**:
observed f₂ = 60.03 (test1) and 48.05 (test5) with published percentile/BCa 90% CIs.
Compare via `--golden-shah-external golden/golden_shah_external.json`.

## Golden values (`golden_polli_slow_similar.json`)

Applies to **polli_slow_similar only** — not the fast HC demo or fast-vs-slow pair.

- 6-unit **polli_slow_similar** reference model fits (Weibull α = 58.97, β = 1.4865, R² = 0.9958,
  AIC = 25.48; etc.) reproduce the **published DDSolver supplementary** (Zhang 2010, Polli 1997).
- Similarity (Ref 6-unit mean vs Test = Unit 6): f2 = 91.02, f1 = 1.60,
  DE% = 69.42 / 70.61, MDT = 13.41 / 13.70 (85% truncation keeps t = 5,10,15,20,25).

## Suggested run

```bash
PYTHONPATH=. python run_benchmark.py \
  --data data \
  --pairs "polli_slow_similar_ref.csv:polli_slow_similar_test.csv,mendyk_reference_N12.csv:mendyk_test_N12.csv,shah_ref.csv:shah_test.csv,shah_ref.csv:shah_borderline_test.csv" \
  --golden golden/golden_polli_slow_similar.json \
  --golden-shah-external golden/golden_shah_external.json
```

The **polli_slow_similar** pair should PASS golden check (tol 0.01).

## Shah 1998 provenance

Canonical source: R package **bootf2** v0.4.1, data object `shah1998`. See `data/shah1998_provenance.json`.
