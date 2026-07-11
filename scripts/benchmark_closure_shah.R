#!/usr/bin/env Rscript
# Shah bootstrap external validation — bootf2 percentile + BCa CIs (est.f2).
suppressPackageStartupMessages({
  if (!requireNamespace("bootf2", quietly = TRUE)) {
    stop("bootf2 package required")
  }
})

args <- commandArgs(trailingOnly = TRUE)
data_dir <- if (length(args) >= 1) args[1] else "data"
n_boots <- if (length(args) >= 2) as.integer(args[2]) else 10000L
seed <- if (length(args) >= 3) as.integer(args[3]) else 42L

run_bootf2_pair <- function(ref_file, test_file, label) {
  ref_df <- read.csv(file.path(data_dir, ref_file), check.names = FALSE)
  test_df <- read.csv(file.path(data_dir, test_file), check.names = FALSE)
  time_col_r <- names(ref_df)[grep("^time", names(ref_df), ignore.case = TRUE)][1]
  time_col_t <- names(test_df)[grep("^time", names(test_df), ignore.case = TRUE)][1]
  names(ref_df)[names(ref_df) == time_col_r] <- "Time"
  names(test_df)[names(test_df) == time_col_t] <- "Time"
  set.seed(seed)
  res <- bootf2::bootf2(
    test_df, ref_df,
    n.boots = n_boots,
    seed = seed,
    regulation = "EMA",
    alpha = 0.05,
    print.report = FALSE,
    output.to.screen = FALSE
  )
  est <- res$boot.summary[res$boot.summary$f2.type == "est.f2", ]
  exp <- res$boot.summary[res$boot.summary$f2.type == "exp.f2", ]
  bc <- res$boot.summary[res$boot.summary$f2.type == "bc.f2", ]
  ci <- res$boot.ci[res$boot.ci$f2.type == "est.f2", ]
  pct <- ci[grepl("Percentile", ci$ci.type, ignore.case = TRUE), ][1, ]
  bca <- ci[grepl("BCa", ci$ci.type, ignore.case = TRUE), ][1, ]
  list(
    label = label,
    f2_observed = est$f2o,
    exp_f2_point = exp$f2o,
    est_boot_mean = est$boot.mean,
    bc_f2_point = bc$f2o,
    percentile_lower = pct$ci.lower,
    percentile_upper = pct$ci.upper,
    bca_lower = bca$ci.lower,
    bca_upper = bca$ci.upper,
    bootf2_version = as.character(packageVersion("bootf2"))
  )
}

cat("=== bootf2 Shah validation (EMA, n=", n_boots, ", seed=", seed, ") ===\n", sep = "")
cat("bootf2 version:", as.character(packageVersion("bootf2")), "\n")
cat("R version:", R.version.string, "\n\n")

pairs <- list(
  run_bootf2_pair("shah_ref.csv", "shah_test.csv", "shah_similar"),
  run_bootf2_pair("shah_ref.csv", "shah_borderline_test.csv", "shah_borderline")
)

for (p in pairs) {
  cat(sprintf("--- %s ---\n", p$label))
  cat(sprintf("  f2_observed=%.4f  exp.f2=%.4f  est.boot.mean=%.4f  bc.f2=%.4f\n",
              p$f2_observed, p$exp_f2_point, p$est_boot_mean, p$bc_f2_point))
  cat(sprintf("  percentile 90%% CI [%.4f, %.4f]\n", p$percentile_lower, p$percentile_upper))
  cat(sprintf("  BCa 90%% CI       [%.4f, %.4f]\n\n", p$bca_lower, p$bca_upper))
}

invisible(pairs)
