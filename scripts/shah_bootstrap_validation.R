#!/usr/bin/env Rscript
# Shah bootstrap validation via bootf2::bootf2 — JSON output for Python orchestrator.
# Usage: Rscript scripts/shah_bootstrap_validation.R [--data data]

args <- commandArgs(trailingOnly = TRUE)
data_dir <- "data"
i <- 1
while (i <= length(args)) {
  if (args[i] == "--data" && i < length(args)) {
    data_dir <- args[i + 1]
    i <- i + 2
  } else {
    i <- i + 1
  }
}

if (!requireNamespace("bootf2", quietly = TRUE)) {
  stop("Install bootf2: install.packages('bootf2')")
}
if (!requireNamespace("jsonlite", quietly = TRUE)) {
  stop("Install jsonlite: install.packages('jsonlite')")
}

run_bootf2_pair <- function(ref_path, test_path, label) {
  ref_df <- read.csv(ref_path, check.names = FALSE)
  test_df <- read.csv(test_path, check.names = FALSE)
  time_col_r <- names(ref_df)[grep("^time", names(ref_df), ignore.case = TRUE)][1]
  time_col_t <- names(test_df)[grep("^time", names(test_df), ignore.case = TRUE)][1]
  names(ref_df)[names(ref_df) == time_col_r] <- "Time"
  names(test_df)[names(test_df) == time_col_t] <- "Time"

  res <- bootf2::bootf2(
    test_df, ref_df,
    n.boots = 10000,
    seed = 42,
    regulation = "EMA",
    alpha = 0.05,
    print.report = FALSE,
    output.to.screen = FALSE
  )

  pick_ci <- function(f2_type, ci_type) {
    rows <- res$boot.ci[res$boot.ci$f2.type == f2_type & grepl(ci_type, res$boot.ci$ci.type, fixed = TRUE), ]
    if (nrow(rows) == 0) return(list(lower = NA_real_, upper = NA_real_))
    list(lower = rows$ci.lower[1], upper = rows$ci.upper[1])
  }

  est <- res$boot.summary[res$boot.summary$f2.type == "est.f2", ]
  bc <- res$boot.summary[res$boot.summary$f2.type == "bc.f2", ]

  pct_est <- pick_ci("est.f2", "Percentile")
  bca_est <- pick_ci("est.f2", "BCa")
  pct_bc <- pick_ci("bc.f2", "Percentile")
  bca_bc <- pick_ci("bc.f2", "BCa")

  list(
    label = label,
    ref_file = basename(ref_path),
    test_file = basename(test_path),
    bootf2_version = as.character(packageVersion("bootf2")),
    est_f2 = list(
      f2_observed = est$f2o[1],
      boot_mean = est$boot.mean[1],
      boot_median = est$boot.median[1],
      percentile_ci_lower = pct_est$lower,
      percentile_ci_upper = pct_est$upper,
      bca_ci_lower = bca_est$lower,
      bca_ci_upper = bca_est$upper
    ),
    bc_f2 = list(
      f2_observed = bc$f2o[1],
      boot_mean = bc$boot.mean[1],
      boot_median = bc$boot.median[1],
      percentile_ci_lower = pct_bc$lower,
      percentile_ci_upper = pct_bc$upper,
      bca_ci_lower = bca_bc$lower,
      bca_ci_upper = bca_bc$upper
    )
  )
}

out <- list(
  generated = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
  bootf2_version = as.character(packageVersion("bootf2")),
  pairs = list(
    shah_similar = run_bootf2_pair(
      file.path(data_dir, "shah_ref.csv"),
      file.path(data_dir, "shah_test.csv"),
      "shah_similar"
    ),
    shah_borderline = run_bootf2_pair(
      file.path(data_dir, "shah_ref.csv"),
      file.path(data_dir, "shah_borderline_test.csv"),
      "shah_borderline"
    )
  )
)

cat(jsonlite::toJSON(out, pretty = TRUE, auto_unbox = TRUE))
