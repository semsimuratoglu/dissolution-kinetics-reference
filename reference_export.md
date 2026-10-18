# Reference engine export

Full export from the independent Python reference engine in this repository.
Generated: 2026-10-18 16:33 UTC

## Canonical dataset registry

| Label | Ref | Test | Role |
|-------|-----|------|------|
| `polli_slow_similar` | `polli_slow_similar_ref.csv` | `polli_slow_similar_test.csv` | Deterministic-equivalence validation vs golden_polli_slow_similar.json |
| `polli_fast_hc_demo` | `polli_fast_hc_demo_ref.csv` | `—` | Hixson-Crowell unconstrained-fit demonstration (kHC≈0.033, R²≈0.83) |
| `polli_fast_vs_slow_dissimilar` | `polli_fast_vs_slow_dissimilar_ref.csv` | `polli_fast_vs_slow_dissimilar_test.csv` | Dissimilarity demonstration pair (f2≈22) |
| `mendyk_high_cv` | `mendyk_reference_N12.csv` | `mendyk_test_N12.csv` | Synthetic high-CV bootstrap anchor |
| `shah_similar` | `shah_ref.csv` | `shah_test.csv` | bootf2::shah1998 test1, f₂≈60.03 |
| `shah_borderline` | `shah_ref.csv` | `shah_borderline_test.csv` | bootf2::shah1998 test5, f₂≈48.05 |

---

## 1. golden_polli_slow_similar.json

```json
{
  "_meta": {
    "polli_set": "polli_slow_similar",
    "description": "Externally-anchored golden for polli_slow_similar only (NOT polli_fast_hc_demo or polli_fast_vs_slow_dissimilar). Model fits on the 6-unit slow metoprolol reference match DDSolver supplementary (Zhang 2010; Polli 1997) to 4 dp. Similarity: Ref(6-unit mean) vs Test(Unit 6). AIC=N*ln(RSS)+2P, FDA/EMA 85% truncation.",
    "data_files": {
      "ref": "polli_slow_similar_ref.csv",
      "test": "polli_slow_similar_test.csv"
    },
    "tolerance_deterministic": 0.01
  },
  "polli_slow_similar": {
    "reference_model_fits": {
      "Weibull": {
        "alpha": 58.9722,
        "beta": 1.4865,
        "R2": 0.9958,
        "AIC": 25.48,
        "RMSE": 2.0976
      },
      "First-order": {
        "k": 0.0668,
        "R2": 0.9405,
        "AIC": 42.35
      },
      "Hixson-Crowell": {
        "kHC": 0.0183,
        "R2": 0.9797,
        "AIC": 34.75
      },
      "Korsmeyer-Peppas_full": {
        "kKP": 12.1171,
        "n": 0.5833,
        "R2": 0.9042,
        "AIC": 47.72
      }
    },
    "similarity_ref_vs_test_unit6": {
      "f2": 91.02,
      "f1": 1.6,
      "DE_ref": 69.42,
      "DE_test": 70.61,
      "MDT_ref": 13.41,
      "MDT_test": 13.7,
      "truncation_points_min": [
        5,
        10,
        15,
        20,
        25
      ]
    }
  }
}
```

## 2. golden_shah_external.json (Noce 2020 / Boddu 2024)

```json
{
  "_meta": {
    "source": "Noce 2020 Table 5; Boddu 2024 Table IV",
    "description": "Externally-anchored Shah dissolution similarity oracle. Observed f2 values are the regulatory mean-profile est.f2 on bootf2::shah1998 (ref vs test1/test5). Bootstrap 90% CIs transcribed from Noce et al. Eur J Pharm Sci 2020 Table 5 (1000 boot) and Boddu et al. AAPS PharmSciTech 2024 Table IV (RLD vs T1/T5). Literature reports percentile and BCa intervals for E(f2)/est.f2; engine does not bootstrap Shah pairs \u2014 CI fields are reference-only (SKIP_NO_BOOT at runtime).",
    "citations": [
      {
        "key": "Noce2020",
        "doi": "10.1016/j.ejps.2020.105259",
        "table": "Table 5",
        "note": "1000 bootstrap replicates"
      },
      {
        "key": "Boddu2024",
        "doi": "10.1208/s12249-023-02710-9",
        "table": "Table IV",
        "note": "RLD vs T1/T5 BCa ranges"
      },
      {
        "key": "Shah1998",
        "doi": "10.1023/A:1011976615750",
        "data": "bootf2::shah1998"
      }
    ],
    "tolerance_f2": 0.01,
    "tolerance_ci": 2.0
  },
  "shah_similar": {
    "pair_label": "shah_similar",
    "ref_file": "shah_ref.csv",
    "test_file": "shah_test.csv",
    "test_arm": "test1",
    "f2_observed": 60.03,
    "f1_observed": 8.73,
    "bootstrap": {
      "percentile_90_ci_lower": 53.2,
      "percentile_90_ci_upper": 67.5,
      "bca_jackknife_ci_lower": 53.9,
      "bca_jackknife_ci_upper": 70.2,
      "boddu_2024_bca_lower_range": [
        53.8,
        55.8
      ],
      "boddu_2024_bca_upper_range": [
        69.3,
        73.0
      ],
      "similar_by_percentile": true
    }
  },
  "shah_borderline": {
    "pair_label": "shah_borderline",
    "ref_file": "shah_ref.csv",
    "test_file": "shah_borderline_test.csv",
    "test_arm": "test5",
    "f2_observed": 48.05,
    "f1_observed": 13.95,
    "bootstrap": {
      "percentile_90_ci_lower": 45.99,
      "percentile_90_ci_upper": 49.96,
      "bca_jackknife_ci_lower": 46.15,
      "bca_jackknife_ci_upper": 50.17,
      "boddu_2024_bca_lower_range": [
        45.9,
        46.3
      ],
      "boddu_2024_bca_upper_range": [
        49.95,
        50.28
      ],
      "similar_by_percentile": false
    }
  }
}
```

## 3. Python reference engine — canonical pairs

Conventions: per-unit fit → mean±SD; RMSE = √(RSS/(n−p)); AIC = N·ln(RSS)+2P;
f₂/f₁ with FDA+EMA 85% truncation; DE%/MDT on full profile. Bootstrap for mendyk only (EMA-primary).

### polli_slow_similar (`polli_slow_similar_ref.csv` vs `polli_slow_similar_test.csv`)

- Time points: 7 | Ref units: 6 | Test units: 1
- Time vector: 5, 10, 15, 20, 25, 30, 45

#### Reference profile — model fits

| Model | Convergence | Parameters (mean±SD) | R² | adj-R² | RMSE | AIC | MSC |
|-------|-------------|----------------------|-----|--------|------|-----|-----|
| Zero-order | 6/6 | k0=2.897461±0.046111 | 0.6119 | 0.5343 | 18.6101 | 55.4534 | 0.6712 |
| First-order | 6/6 | k1=0.066840±0.002813 | 0.9405 | 0.9285 | 7.3154 | 42.3490 | 2.5432 |
| Higuchi | 6/6 | kH=15.828010±0.275807 | 0.8884 | 0.8661 | 10.0192 | 46.7973 | 1.9078 |
| Hixson-Crowell | 6/6 | kHC=0.018343±0.000700 | 0.9797 | 0.9756 | 4.2653 | 34.7482 | 3.6291 |
| Korsmeyer-Peppas (full profile) | 6/6 | kKP=12.117044±1.073377, n=0.583279±0.024442 | 0.9042 | 0.8563 | 10.1576 | 47.7168 | 1.7764 |
| Korsmeyer-Peppas (≤60%) | 2/6 | kKP=3.221528±0.037114, n=1.071515±0.012261 | 0.9984 | n/a | 1.1526 | 4.8519 | 5.1122 |
| | *converged 2/6 units (<=60% window)* | | | | | | |
| Weibull | 6/6 | alpha=58.972229±9.676874, beta=1.486537±0.052062 | 0.9958 | 0.9938 | 2.0976 | 25.4829 | 4.9527 |

#### Test profile — model fits

| Model | Convergence | Parameters (mean±SD) | R² | adj-R² | RMSE | AIC | MSC |
|-------|-------------|----------------------|-----|--------|------|-----|-----|
| Zero-order | 1/1 | k0=2.955814±0.000000 | 0.6443 | 0.5731 | 18.5593 | 55.4359 | 0.7478 |
| First-order | 1/1 | k1=0.068666±0.000000 | 0.9272 | 0.9126 | 8.3974 | 44.3332 | 2.3339 |
| Higuchi | 1/1 | kH=16.118861±0.000000 | 0.8860 | 0.8632 | 10.5073 | 47.4714 | 1.8856 |
| Hixson-Crowell | 1/1 | kHC=0.018781±0.000000 | 0.9731 | 0.9678 | 5.0988 | 37.3483 | 3.3318 |
| Korsmeyer-Peppas (full profile) | 1/1 | kKP=11.807530±0.000000, n=0.595845±0.000000 | 0.9053 | 0.8580 | 10.4883 | 48.1697 | 1.7858 |
| Korsmeyer-Peppas (≤60%) | 0/1 | n/a | n/a | n/a | n/a | n/a | n/a |
| | *converged 0/1 units (<=60% window)* | | | | | | |
| Weibull | 1/1 | alpha=70.259834±0.000000, beta=1.564570±0.000000 | 0.9948 | 0.9922 | 2.4621 | 27.8803 | 4.6843 |

#### Similarity (reference mean vs test mean)

| Metric | FDA | EMA |
|--------|-----|-----|
| f₂ | 91.0165 | 91.0165 |
| f₁ | 1.6033 | 1.6033 |
| DE% (ref) | 69.42 | n/a |
| DE% (test) | 70.61 | n/a |
| MDT (ref) | 13.41 | n/a |
| MDT (test) | 13.70 | n/a |
| Truncation times | 5, 10, 15, 20, 25 | 5, 10, 15, 20, 25 |

#### CV gate / Bootstrap

**CV gate OK — bootstrap not required**

### polli_fast_vs_slow_dissimilar (`polli_fast_vs_slow_dissimilar_ref.csv` vs `polli_fast_vs_slow_dissimilar_test.csv`)

- Time points: 7 | Ref units: 6 | Test units: 6
- Time vector: 5, 10, 15, 20, 25, 30, 45

#### Reference profile — model fits

| Model | Convergence | Parameters (mean±SD) | R² | adj-R² | RMSE | AIC | MSC |
|-------|-------------|----------------------|-----|--------|------|-----|-----|
| Zero-order | 6/6 | k0=3.355407±0.033066 | -2.3610 | -3.0332 | 39.5340 | 66.0207 | -1.4958 |
| First-order | 6/6 | k1=0.153381±0.004381 | 0.9047 | 0.8857 | 6.6594 | 41.0741 | 2.0680 |
| Higuchi | 6/6 | kH=19.250476±0.186555 | 0.2093 | 0.0512 | 19.1639 | 55.8795 | -0.0470 |
| Hixson-Crowell | 6/6 | kHC=0.032749±0.000242 | 0.8294 | 0.7953 | 8.8758 | 45.0794 | 1.4958 |
| Korsmeyer-Peppas (full profile) | 6/6 | kKP=41.536586±1.342065, n=0.260285±0.010008 | 0.6465 | 0.4698 | 14.0513 | 52.2607 | 0.4699 |
| Korsmeyer-Peppas (≤60%) | 0/6 | n/a | n/a | n/a | n/a | n/a | n/a |
| | *converged 0/6 units (<=60% window)* | | | | | | |
| Weibull | 6/6 | alpha=31.341060±2.753291, beta=1.750203±0.048928 | 0.9954 | 0.9931 | 1.5333 | 20.7558 | 4.9706 |

#### Test profile — model fits

| Model | Convergence | Parameters (mean±SD) | R² | adj-R² | RMSE | AIC | MSC |
|-------|-------------|----------------------|-----|--------|------|-----|-----|
| Zero-order | 6/6 | k0=2.897461±0.046111 | 0.6119 | 0.5343 | 18.6101 | 55.4534 | 0.6712 |
| First-order | 6/6 | k1=0.066840±0.002813 | 0.9405 | 0.9285 | 7.3154 | 42.3490 | 2.5432 |
| Higuchi | 6/6 | kH=15.828010±0.275807 | 0.8884 | 0.8661 | 10.0192 | 46.7973 | 1.9078 |
| Hixson-Crowell | 6/6 | kHC=0.018343±0.000700 | 0.9797 | 0.9756 | 4.2653 | 34.7482 | 3.6291 |
| Korsmeyer-Peppas (full profile) | 6/6 | kKP=12.117044±1.073377, n=0.583279±0.024442 | 0.9042 | 0.8563 | 10.1576 | 47.7168 | 1.7764 |
| Korsmeyer-Peppas (≤60%) | 2/6 | kKP=3.221528±0.037114, n=1.071515±0.012261 | 0.9984 | n/a | 1.1526 | 4.8519 | 5.1122 |
| | *converged 2/6 units (<=60% window)* | | | | | | |
| Weibull | 6/6 | alpha=58.972229±9.676874, beta=1.486537±0.052062 | 0.9958 | 0.9938 | 2.0976 | 25.4829 | 4.9527 |

#### Similarity (reference mean vs test mean)

| Metric | FDA | EMA |
|--------|-----|-----|
| f₂ | 26.2179 | 22.4900 |
| f₁ | 32.7581 | 46.2289 |
| DE% (ref) | 85.47 | n/a |
| DE% (test) | 69.42 | n/a |
| MDT (ref) | 6.23 | n/a |
| MDT (test) | 13.41 | n/a |
| Truncation times | 5, 10, 15, 20, 25 | 5, 10, 15 |

#### CV gate / Bootstrap

**CV gate OK — bootstrap not required**

### mendyk_reference_N12_vs_mendyk_test_N12 (`mendyk_reference_N12.csv` vs `mendyk_test_N12.csv`)

- Time points: 6 | Ref units: 12 | Test units: 12
- Time vector: 15, 30, 45, 60, 90, 120

#### Reference profile — model fits

| Model | Convergence | Parameters (mean±SD) | R² | adj-R² | RMSE | AIC | MSC |
|-------|-------------|----------------------|-----|--------|------|-----|-----|
| Zero-order | 12/12 | k0=0.960239±0.063500 | 0.6400 | 0.5500 | 15.7431 | 44.4031 | 0.8296 |
| First-order | 12/12 | k1=0.021244±0.003894 | 0.9790 | 0.9738 | 3.9599 | 28.0938 | 3.5479 |
| Higuchi | 12/12 | kH=8.817844±0.670616 | 0.9353 | 0.9191 | 6.9177 | 34.7052 | 2.4460 |
| Hixson-Crowell | 12/12 | kHC=0.005858±0.000997 | 0.9942 | 0.9927 | 2.0232 | 19.7070 | 4.9457 |
| Korsmeyer-Peppas (full profile) | 12/12 | kKP=6.462301±2.631212, n=0.589167±0.076555 | 0.9599 | 0.9331 | 6.1085 | 34.0224 | 2.5598 |
| Korsmeyer-Peppas (≤60%) | 6/12 | kKP=1.499204±0.415880, n=0.948840±0.059325 | 0.9996 | n/a | 0.4017 | -3.6402 | 7.6102 |
| | *converged 6/12 units (<=60% window)* | | | | | | |
| Weibull | 12/12 | alpha=111.596612±39.316148, beta=1.202724±0.045147 | 0.9950 | 0.9917 | 2.0058 | 19.9585 | 4.9037 |

#### Test profile — model fits

| Model | Convergence | Parameters (mean±SD) | R² | adj-R² | RMSE | AIC | MSC |
|-------|-------------|----------------------|-----|--------|------|-----|-----|
| Zero-order | 12/12 | k0=0.918654±0.067084 | 0.7446 | 0.6808 | 13.4256 | 42.4494 | 1.1716 |
| First-order | 12/12 | k1=0.018930±0.003411 | 0.9776 | 0.9720 | 4.1213 | 28.5733 | 3.4843 |
| Higuchi | 12/12 | kH=8.386073±0.693422 | 0.9260 | 0.9075 | 7.4119 | 35.4557 | 2.3372 |
| Hixson-Crowell | 12/12 | kHC=0.005263±0.000873 | 0.9934 | 0.9918 | 2.1916 | 20.7825 | 4.7828 |
| Korsmeyer-Peppas (full profile) | 12/12 | kKP=5.056001±2.061498, n=0.634983±0.077034 | 0.9651 | 0.9418 | 5.7620 | 33.3236 | 2.6926 |
| Korsmeyer-Peppas (≤60%) | 9/12 | kKP=1.503914±0.809437, n=0.962976±0.116437 | 0.9994 | n/a | 0.4622 | -2.9709 | 7.3290 |
| | *converged 9/12 units (<=60% window)* | | | | | | |
| Weibull | 12/12 | alpha=133.316239±47.114085, beta=1.215994±0.042234 | 0.9958 | 0.9929 | 1.9078 | 19.5601 | 4.9865 |

#### Similarity (reference mean vs test mean)

| Metric | FDA | EMA |
|--------|-----|-----|
| f₂ | 70.9176 | 69.8230 |
| f₁ | 5.6532 | 6.6813 |
| DE% (ref) | 63.52 | n/a |
| DE% (test) | 60.13 | n/a |
| MDT (ref) | 39.85 | n/a |
| MDT (test) | 42.15 | n/a |
| Truncation times | 15, 30, 45, 60, 90, 120 | 15, 30, 45, 60, 90 |

#### CV gate / Bootstrap

**CV gate EXCEEDED — bootstrap EMA (ref: t0 CV=25.73% > 20%; t1 CV=15.33% > 10%; t2 CV=11.95% > 10%; test: t0 CV=28.37% > 20%; t1 CV=16.41% > 10%; t2 CV=11.81% > 10%)**

| Field | Value |
|-------|-------|
| f₂ observed (EMA) | 69.8230 |
| E(f₂) point | 66.1539 |
| E(f₂) median (bootstrap) | 66.1144 |
| bc-f₂ | 73.4921 |
| 90% CI lower | 53.8830 |
| 90% CI upper | 78.9378 |
| similar (CI lower ≥ 50) | True |
| seed | 42 |

### shah_ref_vs_shah_test (`shah_ref.csv` vs `shah_test.csv`)

- Time points: 4 | Ref units: 12 | Test units: 12
- Time vector: 30, 60, 90, 180

#### Reference profile — model fits

| Model | Convergence | Parameters (mean±SD) | R² | adj-R² | RMSE | AIC | MSC |
|-------|-------------|----------------------|-----|--------|------|-----|-----|
| Zero-order | 12/12 | k0=0.641444±0.018595 | 0.2160 | -0.1761 | 22.9048 | 31.4299 | -0.2443 |
| First-order | 12/12 | k1=0.015750±0.001000 | 0.9827 | 0.9740 | 3.3687 | 15.8217 | 3.6578 |
| Higuchi | 12/12 | kH=7.443572±0.218847 | 0.9289 | 0.8933 | 6.8702 | 21.6517 | 2.2003 |
| Hixson-Crowell | 12/12 | kHC=0.004399±0.000269 | 0.9838 | 0.9757 | 2.9439 | 14.4231 | 4.0074 |
| Korsmeyer-Peppas (full profile) | 12/12 | kKP=8.048684±0.875752, n=0.484429±0.022720 | 0.9315 | 0.7946 | 8.2466 | 23.4705 | 1.7456 |
| Korsmeyer-Peppas (≤60%) | 0/12 | n/a | n/a | n/a | n/a | n/a | n/a |
| | *converged 0/12 units (<=60% window)* | | | | | | |
| Weibull | 12/12 | alpha=130.628268±50.057173, beta=1.156082±0.112378 | 0.9945 | 0.9835 | 2.1627 | 12.3024 | 4.5376 |

#### Test profile — model fits

| Model | Convergence | Parameters (mean±SD) | R² | adj-R² | RMSE | AIC | MSC |
|-------|-------------|----------------------|-----|--------|------|-----|-----|
| Zero-order | 12/12 | k0=0.681336±0.024777 | -0.3042 | -0.9563 | 28.1215 | 33.0080 | -0.6924 |
| First-order | 12/12 | k1=0.019429±0.003088 | 0.9750 | 0.9625 | 3.8688 | 16.7895 | 3.3623 |
| Higuchi | 12/12 | kH=7.993334±0.362518 | 0.8394 | 0.7591 | 9.5429 | 24.1803 | 1.5146 |
| Hixson-Crowell | 12/12 | kHC=0.005328±0.000747 | 0.9922 | 0.9883 | 2.1385 | 11.9775 | 4.5653 |
| Korsmeyer-Peppas (full profile) | 12/12 | kKP=11.879931±3.792634, n=0.423574±0.052196 | 0.8877 | 0.6632 | 10.1323 | 25.1902 | 1.2621 |
| Korsmeyer-Peppas (≤60%) | 0/12 | n/a | n/a | n/a | n/a | n/a | n/a |
| | *converged 0/12 units (<=60% window)* | | | | | | |
| Weibull | 12/12 | alpha=135.855112±30.469779, beta=1.233103±0.082723 | 0.9943 | 0.9828 | 2.2124 | 12.3609 | 4.4694 |

#### Similarity (reference mean vs test mean)

| Metric | FDA | EMA |
|--------|-----|-----|
| f₂ | 60.0299 | 57.4692 |
| f₁ | 8.7287 | 11.9811 |
| DE% (ref) | 65.93 | n/a |
| DE% (test) | 71.35 | n/a |
| MDT (ref) | 55.18 | n/a |
| MDT (test) | 48.58 | n/a |
| Truncation times | 30, 60, 90, 180 | 30, 60, 90 |

#### CV gate / Bootstrap

**CV gate OK — bootstrap skipped (Shah deterministic anchor policy)**

### shah_ref_vs_shah_borderline_test (`shah_ref.csv` vs `shah_borderline_test.csv`)

- Time points: 4 | Ref units: 12 | Test units: 12
- Time vector: 30, 60, 90, 180

#### Reference profile — model fits

| Model | Convergence | Parameters (mean±SD) | R² | adj-R² | RMSE | AIC | MSC |
|-------|-------------|----------------------|-----|--------|------|-----|-----|
| Zero-order | 12/12 | k0=0.641444±0.018595 | 0.2160 | -0.1761 | 22.9048 | 31.4299 | -0.2443 |
| First-order | 12/12 | k1=0.015750±0.001000 | 0.9827 | 0.9740 | 3.3687 | 15.8217 | 3.6578 |
| Higuchi | 12/12 | kH=7.443572±0.218847 | 0.9289 | 0.8933 | 6.8702 | 21.6517 | 2.2003 |
| Hixson-Crowell | 12/12 | kHC=0.004399±0.000269 | 0.9838 | 0.9757 | 2.9439 | 14.4231 | 4.0074 |
| Korsmeyer-Peppas (full profile) | 12/12 | kKP=8.048684±0.875752, n=0.484429±0.022720 | 0.9315 | 0.7946 | 8.2466 | 23.4705 | 1.7456 |
| Korsmeyer-Peppas (≤60%) | 0/12 | n/a | n/a | n/a | n/a | n/a | n/a |
| | *converged 0/12 units (<=60% window)* | | | | | | |
| Weibull | 12/12 | alpha=130.628268±50.057173, beta=1.156082±0.112378 | 0.9945 | 0.9835 | 2.1627 | 12.3024 | 4.5376 |

#### Test profile — model fits

| Model | Convergence | Parameters (mean±SD) | R² | adj-R² | RMSE | AIC | MSC |
|-------|-------------|----------------------|-----|--------|------|-----|-----|
| Zero-order | 12/12 | k0=0.699822±0.012975 | -0.8524 | -1.7787 | 32.1188 | 34.1432 | -1.1071 |
| First-order | 12/12 | k1=0.022041±0.001016 | 0.9676 | 0.9514 | 4.2472 | 17.8199 | 2.9737 |
| Higuchi | 12/12 | kH=8.286116±0.152971 | 0.7486 | 0.6229 | 11.7858 | 26.0919 | 0.9057 |
| Hixson-Crowell | 12/12 | kHC=0.005992±0.000248 | 0.9720 | 0.9580 | 3.9024 | 17.2175 | 3.1243 |
| Korsmeyer-Peppas (full profile) | 12/12 | kKP=15.465939±1.404724, n=0.367023±0.018830 | 0.8531 | 0.5593 | 11.0870 | 25.9743 | 0.9351 |
| Korsmeyer-Peppas (≤60%) | 0/12 | n/a | n/a | n/a | n/a | n/a | n/a |
| | *converged 0/12 units (<=60% window)* | | | | | | |
| Weibull | 12/12 | alpha=110.809981±50.076544, beta=1.208170±0.109067 | 0.9838 | 0.9514 | 3.5424 | 16.3199 | 3.3487 |

#### Similarity (reference mean vs test mean)

| Metric | FDA | EMA |
|--------|-----|-----|
| f₂ | 48.0521 | 45.2334 |
| f₁ | 13.9500 | 19.5711 |
| DE% (ref) | 65.93 | n/a |
| DE% (test) | 73.64 | n/a |
| MDT (ref) | 55.18 | n/a |
| MDT (test) | 45.52 | n/a |
| Truncation times | 30, 60, 90, 180 | 30, 60, 90 |

#### CV gate / Bootstrap

**CV gate OK — bootstrap skipped (Shah deterministic anchor policy)**

---

## 4. Golden / external validation

**polli_slow_similar GoldenCheck:** PASS

| Quantity | Engine | Golden | Abs diff | Status |
|----------|--------|--------|----------|--------|
| polli_slow_similar.ref.Weibull.alpha | 58.9722 | 58.9722 | 2.915e-05 | PASS |
| polli_slow_similar.ref.Weibull.beta | 1.4865 | 1.4865 | 3.678e-05 | PASS |
| polli_slow_similar.ref.Weibull.R2 | 0.9958 | 0.9958 | 4.627e-05 | PASS |
| polli_slow_similar.ref.Weibull.AIC | 25.4829 | 25.4800 | 0.002934 | PASS |
| polli_slow_similar.ref.Weibull.RMSE | 2.0976 | 2.0976 | 5.262e-07 | PASS |
| polli_slow_similar.ref.First-order.k | 0.0668 | 0.0668 | 3.96e-05 | PASS |
| polli_slow_similar.ref.First-order.R2 | 0.9405 | 0.9405 | 4.268e-05 | PASS |
| polli_slow_similar.ref.First-order.AIC | 42.3490 | 42.3500 | 0.0009955 | PASS |
| polli_slow_similar.ref.Hixson-Crowell.kHC | 0.0183 | 0.0183 | 4.322e-05 | PASS |
| polli_slow_similar.ref.Hixson-Crowell.R2 | 0.9797 | 0.9797 | 1.637e-05 | PASS |
| polli_slow_similar.ref.Hixson-Crowell.AIC | 34.7482 | 34.7500 | 0.001807 | PASS |
| polli_slow_similar.ref.Korsmeyer-Peppas_full.kKP | 12.1170 | 12.1171 | 5.57e-05 | PASS |
| polli_slow_similar.ref.Korsmeyer-Peppas_full.n | 0.5833 | 0.5833 | 2.127e-05 | PASS |
| polli_slow_similar.ref.Korsmeyer-Peppas_full.R2 | 0.9042 | 0.9042 | 2.473e-05 | PASS |
| polli_slow_similar.ref.Korsmeyer-Peppas_full.AIC | 47.7168 | 47.7200 | 0.003232 | PASS |
| polli_slow_similar.similarity.f2 | 91.0165 | 91.0200 | 0.0035 | PASS |
| polli_slow_similar.similarity.f1 | 1.6033 | 1.6000 | 0.003336 | PASS |
| polli_slow_similar.similarity.DE_ref | 69.4222 | 69.4200 | 0.002222 | PASS |
| polli_slow_similar.similarity.DE_test | 70.6056 | 70.6100 | 0.004444 | PASS |
| polli_slow_similar.similarity.MDT_ref | 13.4125 | 13.4100 | 0.002538 | PASS |
| polli_slow_similar.similarity.MDT_test | 13.6970 | 13.7000 | 0.002956 | PASS |
| polli_slow_similar.similarity.truncation_points_min | [5.0, 10.0, 15.0, 20.0, 25.0] | [5, 10, 15, 20, 25] | N/A | PASS |

**ShahExternalCheck (literature):** PASS

| Quantity | Engine | Golden | Abs diff | Status |
|----------|--------|--------|----------|--------|
| shah_similar.f2_observed | 60.0299 | 60.0300 | 0.0001233 | PASS |
| shah_similar.f1_observed | 8.7287 | 8.7300 | 0.001278 | PASS |
| shah_similar.boot.percentile_90_ci_lower | None | 53.2000 | N/A | SKIP_NO_BOOT |
| shah_similar.boot.percentile_90_ci_upper | None | 67.5000 | N/A | SKIP_NO_BOOT |
| shah_borderline.f2_observed | 48.0521 | 48.0500 | 0.002052 | PASS |
| shah_borderline.f1_observed | 13.9500 | 13.9500 | 4.496e-05 | PASS |
| shah_borderline.boot.percentile_90_ci_lower | None | 45.9900 | N/A | SKIP_NO_BOOT |
| shah_borderline.boot.percentile_90_ci_upper | None | 49.9600 | N/A | SKIP_NO_BOOT |

## 5. Shah borderline summary

- **f₂ (FDA/EMA):** 48.0521 (literature ≈48.05)
- **CV gate:** CV gate OK — bootstrap skipped (Shah deterministic anchor policy)
- **Bootstrap triggered:** False
- **Verdict:** borderline similarity (point estimate near 50; literature percentile CI upper < 50 for test5)

## 6. R cross-check

Run:

```bash
Rscript reference.R --data data \
  --pairs "polli_slow_similar_ref.csv:polli_slow_similar_test.csv,mendyk_reference_N12.csv:mendyk_test_N12.csv,shah_ref.csv:shah_test.csv,shah_ref.csv:shah_borderline_test.csv" \
  --validate results/python_reference_export.json \
  --golden golden/golden_polli_slow_similar.json
```

Shah export: `Rscript scripts/export_shah1998.R --out-dir data`

