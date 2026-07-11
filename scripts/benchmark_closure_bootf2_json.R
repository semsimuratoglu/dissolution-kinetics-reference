#!/usr/bin/env Rscript
# JSON bootf2 est.f2 CIs for benchmark_closure.py (stdout only).
suppressPackageStartupMessages({
  library(bootf2)
  library(jsonlite)
})

args <- commandArgs(trailingOnly = TRUE)
data_dir <- if (length(args) >= 1) args[1] else "data"
n_boots <- if (length(args) >= 2) as.integer(args[2]) else 10000L
seed <- if (length(args) >= 3) as.integer(args[3]) else 42L

run_pair <- function(ref_file, test_file, label) {
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
  ci <- res$boot.ci[res$boot.ci$f2.type == "est.f2", ]
  pct <- ci[grepl("Percentile", ci$ci.type, ignore.case = TRUE), ][1, ]
  bca <- ci[grepl("BCa", ci$ci.type, ignore.case = TRUE), ][1, ]
  list(
    label = label,
    f2_obs = est$f2o,
    exp_f2 = exp$f2o,
    pct_l = pct$ci.lower,
    pct_u = pct$ci.upper,
    bca_l = bca$ci.lower,
    bca_u = bca$ci.upper
  )
}

out <- list(
  run_pair("shah_ref.csv", "shah_test.csv", "shah_similar"),
  run_pair("shah_ref.csv", "shah_borderline_test.csv", "shah_borderline")
)
cat(toJSON(out, auto_unbox = TRUE))
