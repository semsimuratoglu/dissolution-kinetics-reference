#!/usr/bin/env Rscript
# Export full R-engine reference results as JSON (stdout).
# Independent minpack.lm fits + unit-level EMA-primary bootstrap (seed=42, n=10000).

suppressPackageStartupMessages({
  library(minpack.lm)
  library(jsonlite)
})

args <- commandArgs(trailingOnly = TRUE)
data_dir <- if (length(args) >= 1) args[1] else "data"
ref_r_path <- if (length(args) >= 2) args[2] else {
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  script_dir <- if (length(file_arg)) {
    dirname(sub("^--file=", "", file_arg[1]))
  } else {
    "."
  }
  file.path(script_dir, "..", "reference.R")
}

BOOT_ITER <- 10000L
BOOT_SEED <- 42L

source_reference <- function() {
  ref_path <- ref_r_path
  if (!file.exists(ref_path)) {
    stop("reference.R not found at ", ref_path)
  }
  lines <- readLines(ref_path, warn = FALSE)
  lines <- lines[!grepl("^main\\(\\)", lines)]
  tf <- tempfile(fileext = ".R")
  writeLines(lines, tf)
  source(tf, local = FALSE)
  unlink(tf)
}

# Re-use reference.R implementations without running main().
tryCatch(
  source_reference(),
  error = function(e) stop("Failed to source reference.R: ", conditionMessage(e))
)

DATASETS <- list(
  list(
    label = "polli_slow_similar",
    ref = "polli_slow_similar_ref.csv",
    test = "polli_slow_similar_test.csv",
    time_unit = "min"
  ),
  list(
    label = "polli_fast_vs_slow_dissimilar",
    ref = "polli_fast_vs_slow_dissimilar_ref.csv",
    test = "polli_fast_vs_slow_dissimilar_test.csv",
    time_unit = "min"
  ),
  list(
    label = "mendyk",
    ref = "mendyk_reference_N12.csv",
    test = "mendyk_test_N12.csv",
    time_unit = "min"
  ),
  list(
    label = "shah_similar",
    ref = "shah_ref.csv",
    test = "shah_test.csv",
    time_unit = "min"
  ),
  list(
    label = "shah_borderline",
    ref = "shah_ref.csv",
    test = "shah_borderline_test.csv",
    time_unit = "min"
  ),
  list(
    label = "tsong",
    ref = "tsong_ref.csv",
    test = "tsong_test.csv",
    time_unit = "min"
  ),
  list(
    label = "chow",
    ref = "chow_ref.csv",
    test = "chow_test.csv",
    time_unit = "hours"
  )
)

MODEL_KEYS <- names(MODEL_REGISTRY)

first_point_cv <- function(ref_units, test_units) {
  cv_t0 <- function(units) {
    if (nrow(units) < 2) return(NA_real_)
    means <- colMeans(units)
    sds <- apply(units, 2, sd)
    if (means[1] == 0) return(NA_real_)
    100 * sds[1] / abs(means[1])
  }
  ref_cv <- cv_t0(ref_units)
  test_cv <- cv_t0(test_units)
  if (is.na(ref_cv)) return(test_cv)
  if (is.na(test_cv)) return(ref_cv)
  max(ref_cv, test_cv)
}

serialize_model <- function(summary) {
  params <- summary$param_mean
  if (length(params) == 0) params <- list()
  out <- list(
    converged = paste0(summary$n_converged, "/", summary$n_total),
    n_converged = summary$n_converged,
    n_total = summary$n_total,
    params = params,
    gof = list(
      R2 = summary$gof_mean$R2,
      adjR2 = summary$gof_mean$adj_R2,
      RMSE = summary$gof_mean$RMSE,
      AIC = summary$gof_mean$AIC
    )
  )
  out
}

best_model_name <- function(summaries) {
  best_key <- NA_character_
  best_aic <- Inf
  for (mk in MODEL_KEYS) {
    if (mk == "korsmeyer_peppas_60") next
    sm <- summaries[[mk]]
    if (sm$n_converged <= 0.5 * sm$n_total) next
    aic_val <- sm$gof_mean$AIC
    if (is.null(aic_val) || is.na(aic_val)) next
    if (aic_val < best_aic) {
      best_aic <- aic_val
      best_key <- mk
    }
  }
  best_key
}

export_dataset <- function(spec) {
  ref_path <- file.path(data_dir, spec$ref)
  test_path <- file.path(data_dir, spec$test)
  ref <- load_csv(ref_path)
  test <- load_csv(test_path)

  ref_models <- lapply(MODEL_KEYS, function(mk) {
    fit_model_all_units(ref$times, ref$units, mk)
  })
  names(ref_models) <- MODEL_KEYS

  test_models <- lapply(MODEL_KEYS, function(mk) {
    fit_model_all_units(test$times, test$units, mk)
  })
  names(test_models) <- MODEL_KEYS

  metrics <- compute_profile_metrics(ref$times, ref$units, test$units)
  cv_ref <- cv_gate(ref$units)
  cv_test <- cv_gate(test$units)
  high_cv <- cv_ref$exceeded || cv_test$exceeded

  boot <- NULL
  if (high_cv && grepl("mendyk", spec$label, fixed = TRUE)) {
    boot <- python_style_bootstrap(ref$times, ref$units, test$units)
  }

  list(
    label = spec$label,
    ref_file = spec$ref,
    test_file = spec$test,
    time_unit = spec$time_unit,
    n_timepoints = length(ref$times),
    n_ref_units = nrow(ref$units),
    n_test_units = nrow(test$units),
    similarity = list(
      f2_primary = metrics$f2_primary,
      f1_primary = metrics$f1_primary,
      primary_rule = metrics$primary,
      first_point_CV = first_point_cv(ref$units, test$units),
      high_variability = high_cv
    ),
    exposure = list(
      DE_ref = metrics$DE_ref,
      DE_test = metrics$DE_test,
      MDT_ref = metrics$MDT_ref,
      MDT_test = metrics$MDT_test
    ),
    bootstrap = if (!is.null(boot)) {
      list(
        E_f2 = boot$E_f2,
        ci_lower = boot$ci_lower,
        ci_upper = boot$ci_upper,
        decision = if (boot$similar) "similar" else "not similar"
      )
    } else {
      NULL
    },
    ref_models = lapply(ref_models, serialize_model),
    test_models = lapply(test_models, serialize_model),
    best_fit_ref = best_model_name(ref_models)
  )
}

versions <- list(
  python = NA_character_,
  numpy = NA_character_,
  scipy = NA_character_,
  R = paste(R.version$major, R.version$minor, sep = "."),
  minpack_lm = as.character(packageVersion("minpack.lm")),
  bootf2 = tryCatch(as.character(packageVersion("bootf2")), error = function(e) NA_character_),
  export_date = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z", tz = "UTC")
)

out <- list(
  versions = versions,
  datasets = lapply(DATASETS, export_dataset)
)

cat(toJSON(out, auto_unbox = TRUE, digits = 10, na = "null"))
