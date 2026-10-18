# Platform outputs (for independent reproduction)

Output files exported by the ReleaseKinetix platform after the EMA 85% truncation fix.
No platform code is included. The files are byte-identical to the platform exports; only
the file names were changed. They are read by `scripts/compare/build_equivalence_tables.py`
and `scripts/compare/shah_bootstrap_crosscheck.py`.

| File | Original file name | Platform commit | Export date | SHA-256 |
|------|--------------------|-----------------|-------------|---------|
| `platform_results.json` | `results_releasekinetix.json` | `11553fd79ac4e3741444e682719ebaadf2adeae9` | 2026-10-18 | `24122a112272473fb6b2cb10c5ba30034d8d84c66c99539447507dc3aef16341` |
| `platform_shah_bootstrap.json` | `shah_bootstrap_site.json` | `69054b88f76315d08a61f9af60c2e3fe16036579` | 2026-10-18 (`generated_at` 2026-10-18T16:17:50.947Z) | `d545290d54972afd1984caff023a5ef7754ee8d1656dc6cc95b31421dee8775b` |

Platform engine code: tag math-core-v1.1.0 (commit f1e6c8a7822ba99c295640feec671d517c28fe5f); engine code is identical at the export commits listed above.

Contents:

- `platform_results.json`: seven benchmark pairs; reference-arm model fits, AIC candidate sets,
  f2 (FDA, EMA, primary), f1, DE%, MDT, first-point CV, and the mendyk bootstrap block.
- `platform_shah_bootstrap.json`: shah_similar and shah_borderline; unit-level bootstrap of E(f2)
  (10,000 iterations, 90% percentile CI, dynamic EMA truncation per draw) at seed 42 and for
  seeds 42–66.

Verify: `shasum -a 256 platform_outputs/*.json`
