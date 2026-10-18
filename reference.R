#!/usr/bin/env Rscript
# Independent R cross-check engine for dissolution analytics benchmarking.
# Mirrors the Python refeng module (models, GOF, metrics, bootstrap).
#
# Usage:
#   Rscript reference.R --data data \
#     --pairs "polli_slow_similar_ref.csv:polli_slow_similar_test.csv,mendyk_reference_N12.csv:mendyk_test_N12.csv,shah_ref.csv:shah_test.csv"
#
#   Rscript reference.R --data data --pairs ... --validate results/python_reference_export.json \
#     --golden golden/golden_polli_slow_similar.json
#
# Dependencies: minpack.lm, bootf2 (see renv-versions.txt)

args <- commandArgs(trailingOnly = TRUE)

BOOT_ITER <- 10000L
BOOT_SEED <- 42L

parse_args <- function(args) {
  data_dir <- "data"
  pairs <- paste(
    "polli_slow_similar_ref.csv:polli_slow_similar_test.csv",
    "mendyk_reference_N12.csv:mendyk_test_N12.csv",
    "shah_ref.csv:shah_test.csv",
    sep = ","
  )
  validate_json <- NULL
  golden_json <- NULL
  i <- 1
  while (i <= length(args)) {
    if (args[i] == "--data" && i < length(args)) {
      data_dir <- args[i + 1]
      i <- i + 2
    } else if (args[i] == "--pairs" && i < length(args)) {
      pairs <- args[i + 1]
      i <- i + 2
    } else if (args[i] == "--validate" && i < length(args)) {
      validate_json <- args[i + 1]
      i <- i + 2
    } else if (args[i] == "--golden" && i < length(args)) {
      golden_json <- args[i + 1]
      i <- i + 2
    } else {
      i <- i + 1
    }
  }
  list(
    data_dir = data_dir,
    pairs = strsplit(pairs, ",")[[1]],
    validate_json = validate_json,
    golden_json = golden_json
  )
}

load_csv <- function(path) {
  df <- read.csv(path, check.names = FALSE)
  time_col <- names(df)[grep("^time", names(df), ignore.case = TRUE)][1]
  if (is.na(time_col)) stop("Expected a Time column in ", path)
  rep_cols <- setdiff(names(df), time_col)
  if (length(rep_cols) == 0) stop("Expected replicate columns in ", path)
  times <- as.numeric(df[[time_col]])
  units <- t(as.matrix(df[, rep_cols, drop = FALSE]))
  ord <- order(times)
  list(times = times[ord], units = units[, ord, drop = FALSE])
}

cv_gate <- function(units) {
  if (nrow(units) < 2) {
    return(list(exceeded = FALSE, message = ""))
  }
  means <- colMeans(units)
  sds <- apply(units, 2, sd)
  cv <- 100 * sds / abs(means)
  cv[is.nan(cv)] <- 0
  flags <- c()
  if (length(cv) >= 1 && cv[1] > 20) {
    flags <- c(flags, sprintf("t0 CV=%.2f%% > 20%%", cv[1]))
  }
  if (length(cv) >= 2) {
    for (j in 2:length(cv)) {
      if (cv[j] > 10) flags <- c(flags, sprintf("t%d CV=%.2f%% > 10%%", j - 1, cv[j]))
    }
  }
  list(exceeded = length(flags) > 0, message = paste(flags, collapse = "; "))
}

gof_metrics <- function(observed, predicted, n_params) {
  rss <- sum((observed - predicted)^2)
  tss <- sum((observed - mean(observed))^2)
  n <- length(observed)
  r2 <- if (tss == 0) {
    if (rss == 0) 1 else 0
  } else {
    1 - rss / tss
  }
  adj_r2 <- if (n - n_params - 1 > 0) {
    1 - (1 - r2) * (n - 1) / (n - n_params - 1)
  } else {
    NA_real_
  }
  df <- n - n_params
  rmse <- if (df > 0) sqrt(rss / df) else NA_real_
  rss_floor <- max(rss, .Machine$double.xmin)
  tss_floor <- max(tss, .Machine$double.xmin)
  aic <- n * log(rss_floor) + 2 * n_params
  msc <- log(tss_floor / rss_floor) - 2 * n_params / n
  list(R2 = r2, adj_R2 = adj_r2, RMSE = rmse, AIC = aic, MSC = msc, RSS = rss)
}

fit_zero_order_unit <- function(times, y) {
  dat <- data.frame(t = times, Q = y)
  fit <- tryCatch(lm(Q ~ 0 + t, data = dat), error = function(e) NULL)
  if (is.null(fit)) return(NULL)
  k0 <- unname(coef(fit)[1])
  pred <- k0 * times
  list(obs = y, pred = pred, params = list(k0 = k0), n_params = 1L)
}

fit_higuchi_unit <- function(times, y) {
  st <- sqrt(pmax(times, 0))
  dat <- data.frame(st = st, Q = y)
  fit <- tryCatch(lm(Q ~ 0 + st, data = dat), error = function(e) NULL)
  if (is.null(fit)) return(NULL)
  kH <- unname(coef(fit)[1])
  pred <- kH * st
  list(obs = y, pred = pred, params = list(kH = kH), n_params = 1L)
}

fit_first_order_unit <- function(times, y) {
  dat <- data.frame(t = times, Q = y)
  fit <- tryCatch(
    minpack.lm::nlsLM(
      Q ~ 100 * (1 - exp(-k1 * t)),
      data = dat,
      start = list(k1 = 0.05),
      lower = c(k1 = 1e-8)
    ),
    error = function(e) NULL
  )
  if (is.null(fit)) return(NULL)
  pred <- as.numeric(predict(fit, newdata = dat))
  list(obs = y, pred = pred, params = as.list(coef(fit)), n_params = 1L)
}

fit_hixson_crowell_unit <- function(times, y) {
  dat <- data.frame(t = times, Q = y)
  fit <- tryCatch(
    minpack.lm::nlsLM(
      Q ~ 100 * (1 - (1 - kHC * t)^3),
      data = dat,
      start = list(kHC = 0.01),
      lower = c(kHC = 1e-8),
      upper = c(kHC = 1.0)
    ),
    error = function(e) NULL
  )
  if (is.null(fit)) return(NULL)
  pred <- as.numeric(predict(fit, newdata = dat))
  list(obs = y, pred = pred, params = as.list(coef(fit)), n_params = 1L)
}

fit_kp_full_unit <- function(times, y) {
  dat <- data.frame(t = times, Q = y)
  fit <- tryCatch(
    minpack.lm::nlsLM(
      Q ~ kKP * t^n,
      data = dat,
      start = list(kKP = 10, n = 0.5),
      lower = c(kKP = 1e-8, n = 0),
      upper = c(kKP = Inf, n = 10)
    ),
    error = function(e) NULL
  )
  if (is.null(fit)) return(NULL)
  pred <- as.numeric(predict(fit, newdata = dat))
  list(obs = y, pred = pred, params = as.list(coef(fit)), n_params = 2L)
}

fit_kp60_unit <- function(times, y) {
  keep <- y <= 60
  if (sum(keep) < 3) return(NULL)
  t_sub <- times[keep]
  y_sub <- y[keep]
  dat <- data.frame(t = t_sub, Q = y_sub)
  fit <- tryCatch(
    minpack.lm::nlsLM(
      Q ~ kKP * t^n,
      data = dat,
      start = list(kKP = 10, n = 0.5),
      lower = c(kKP = 1e-8, n = 0),
      upper = c(kKP = Inf, n = 10)
    ),
    error = function(e) NULL
  )
  if (is.null(fit)) return(NULL)
  pred <- as.numeric(predict(fit, newdata = dat))
  list(obs = y_sub, pred = pred, params = as.list(coef(fit)), n_params = 2L)
}

fit_weibull_unit <- function(times, y) {
  dat <- data.frame(t = times, Q = y)
  # Log–log linearization seeds (match Python/TS fitWeibullNLS)
  alpha_seed <- 50
  beta_seed <- 1.5
  ok <- is.finite(dat$t) & is.finite(dat$Q) & dat$t > 0 & dat$Q > 0 & dat$Q < 100
  if (sum(ok) >= 2) {
    lx <- log(dat$t[ok])
    ly <- log(-log(1 - dat$Q[ok] / 100))
    if (all(is.finite(lx)) && all(is.finite(ly))) {
      lin <- coef(lm(ly ~ lx))
      beta_seed <- min(max(as.numeric(lin[2]), 0.05), 10)
      k_seed <- exp(as.numeric(lin[1]))
      if (is.finite(k_seed) && k_seed > 0) alpha_seed <- min(max(1 / k_seed, 1e-6), 1e9)
    }
  }
  fit <- tryCatch(
    minpack.lm::nlsLM(
      Q ~ 100 * (1 - exp(-(pmax(t, 0)^beta) / alpha)),
      data = dat,
      start = list(alpha = alpha_seed, beta = beta_seed),
      lower = c(alpha = 1e-8, beta = 1e-8),
      upper = c(alpha = Inf, beta = 10)
    ),
    error = function(e) NULL
  )
  if (is.null(fit)) return(NULL)
  pred <- as.numeric(predict(fit, newdata = dat))
  list(obs = y, pred = pred, params = as.list(coef(fit)), n_params = 2L)
}

MODEL_REGISTRY <- list(
  zero_order = list(label = "Zero-order", fit_fn = fit_zero_order_unit, n_params = 1L),
  first_order = list(label = "First-order", fit_fn = fit_first_order_unit, n_params = 1L),
  higuchi = list(label = "Higuchi", fit_fn = fit_higuchi_unit, n_params = 1L),
  hixson_crowell = list(label = "Hixson-Crowell", fit_fn = fit_hixson_crowell_unit, n_params = 1L),
  korsmeyer_peppas_full = list(label = "K-P full", fit_fn = fit_kp_full_unit, n_params = 2L),
  korsmeyer_peppas_60 = list(label = "K-P <=60%", fit_fn = fit_kp60_unit, n_params = 2L),
  weibull = list(label = "Weibull", fit_fn = fit_weibull_unit, n_params = 2L)
)

fit_model_all_units <- function(times, units, model_key) {
  spec <- MODEL_REGISTRY[[model_key]]
  n_units <- nrow(units)
  params_list <- list()
  gof_list <- list()
  for (i in seq_len(n_units)) {
    res <- spec$fit_fn(times, units[i, ])
    if (is.null(res)) next
    gof <- gof_metrics(res$obs, res$pred, spec$n_params)
    params_list[[length(params_list) + 1]] <- res$params
    gof_list[[length(gof_list) + 1]] <- gof
  }
  n_conv <- length(params_list)
  if (n_conv == 0) {
    return(list(
      model_key = model_key,
      label = spec$label,
      n_converged = 0L,
      n_total = n_units,
      param_mean = list(),
      param_sd = list(),
      gof_mean = list()
    ))
  }
  param_mat <- do.call(rbind, lapply(params_list, function(x) unlist(x)))
  if (is.vector(param_mat)) param_mat <- matrix(param_mat, nrow = 1)
  param_mean <- colMeans(param_mat)
  param_sd <- if (n_conv > 1) apply(param_mat, 2, sd) else rep(0, ncol(param_mat))
  gof_mat <- do.call(rbind, lapply(gof_list, function(x) unlist(x)))
  gof_mean <- colMeans(gof_mat, na.rm = TRUE)
  list(
    model_key = model_key,
    label = spec$label,
    n_converged = n_conv,
    n_total = n_units,
    param_mean = as.list(param_mean),
    param_sd = as.list(param_sd),
    gof_mean = as.list(gof_mean)
  )
}

fda_truncate <- function(times, ref_mean, test_mean) {
  both <- (ref_mean > 85) & (test_mean > 85)
  if (any(both)) {
    last <- which(both)[1]
    idx <- seq_len(last)
  } else {
    idx <- seq_along(times)
  }
  list(times = times[idx], ref = ref_mean[idx], test = test_mean[idx], idx = idx)
}

ema_truncate <- function(times, ref_mean, test_mean) {
  any_above <- (ref_mean > 85) | (test_mean > 85)
  if (any(any_above)) {
    last <- which(any_above)[1]
    idx <- seq_len(last)
  } else {
    idx <- seq_along(times)
  }
  list(times = times[idx], ref = ref_mean[idx], test = test_mean[idx], idx = idx)
}

f2 <- function(ref, test) {
  if (length(ref) == 0) return(NA_real_)
  msq <- mean((ref - test)^2)
  50 * log10(100 / sqrt(1 + msq))
}

f1 <- function(ref, test) {
  ref_sum <- sum(ref)
  if (ref_sum == 0) return(NA_real_)
  100 * sum(abs(ref - test)) / ref_sum
}

de_percent <- function(times, means) {
  if (length(times) == 0) return(NA_real_)
  t_ext <- c(0, times)
  q_ext <- c(0, means)
  auc <- sum(diff(t_ext) * (head(q_ext, -1) + tail(q_ext, -1)) / 2)
  t_last <- times[length(times)]
  if (t_last == 0) return(NA_real_)
  auc / (100 * t_last) * 100
}

mdt <- function(times, means) {
  if (length(times) == 0) return(NA_real_)
  t_ext <- c(0, times)
  q_ext <- c(0, means)
  dQ <- diff(q_ext)
  t_mid <- (head(t_ext, -1) + tail(t_ext, -1)) / 2
  dQ_sum <- sum(dQ)
  if (dQ_sum == 0) return(NA_real_)
  sum(t_mid * dQ) / dQ_sum
}

compute_profile_metrics <- function(times, ref_units, test_units) {
  ref_mean <- colMeans(ref_units)
  test_mean <- colMeans(test_units)
  fda <- fda_truncate(times, ref_mean, test_mean)
  ema <- ema_truncate(times, ref_mean, test_mean)
  slices_diverge <- !identical(fda$idx, ema$idx)
  primary <- if (slices_diverge) "EMA" else "FDA"
  list(
    f2_fda = f2(fda$ref, fda$test),
    f1_fda = f1(fda$ref, fda$test),
    f2_ema = f2(ema$ref, ema$test),
    f1_ema = f1(ema$ref, ema$test),
    primary = primary,
    f2_primary = if (primary == "EMA") f2(ema$ref, ema$test) else f2(fda$ref, fda$test),
    f1_primary = if (primary == "EMA") f1(ema$ref, ema$test) else f1(fda$ref, fda$test),
    DE_ref = de_percent(times, ref_mean),
    DE_test = de_percent(times, test_mean),
    MDT_ref = mdt(times, ref_mean),
    MDT_test = mdt(times, test_mean),
    truncated_times_fda = fda$times,
    truncated_times_ema = ema$times
  )
}

run_bootf2_ema <- function(ref_path, test_path, n_boots = BOOT_ITER, seed = BOOT_SEED) {
  ref_df <- read.csv(ref_path, check.names = FALSE)
  test_df <- read.csv(test_path, check.names = FALSE)
  time_col_r <- names(ref_df)[grep("^time", names(ref_df), ignore.case = TRUE)][1]
  time_col_t <- names(test_df)[grep("^time", names(test_df), ignore.case = TRUE)][1]
  names(ref_df)[names(ref_df) == time_col_r] <- "Time"
  names(test_df)[names(test_df) == time_col_t] <- "Time"
  res <- bootf2::bootf2(
    test_df, ref_df,
    n.boots = n_boots,
    seed = seed,
    regulation = "EMA",
    alpha = 0.05,
    print.report = FALSE,
    output.to.screen = FALSE
  )
  est_row <- res$boot.summary[res$boot.summary$f2.type == "est.f2", ]
  bc_row <- res$boot.summary[res$boot.summary$f2.type == "bc.f2", ]
  ci_row <- res$boot.ci[res$boot.ci$f2.type == "bc.f2" & res$boot.ci$ci.type == "Percentile (Type 7)", ]
  if (nrow(ci_row) == 0) {
    ci_row <- res$boot.ci[res$boot.ci$f2.type == "bc.f2" & grepl("Percentile", res$boot.ci$ci.type), ][1, ]
  }
  list(
    est_f2_boot_mean = est_row$boot.mean,
    est_f2_boot_median = est_row$boot.median,
    bc_f2 = bc_row$f2o,
    bc_f2_ci_lower = ci_row$ci.lower,
    bc_f2_ci_upper = ci_row$ci.upper,
    bc_f2_similar = ci_row$ci.lower >= 50,
    f2_observed = est_row$f2o
  )
}

expected_f2_analytical <- function(times, ref_units, test_units) {
  ref_mean <- colMeans(ref_units)
  test_mean <- colMeans(test_units)
  fda <- fda_truncate(times, ref_mean, test_mean)
  ema <- ema_truncate(times, ref_mean, test_mean)
  slices_diverge <- !identical(fda$idx, ema$idx)
  kept <- if (slices_diverge) ema$idx else fda$idx
  ref_t <- if (slices_diverge) ema$ref else fda$ref
  test_t <- if (slices_diverge) ema$test else fda$test
  n_ref <- nrow(ref_units)
  n_test <- nrow(test_units)
  P <- length(kept)
  if (P < 3) return(NA_real_)
  diff_sq <- sum((ref_t - test_t)^2)
  var_sum <- 0
  for (j in kept) {
    if (n_ref > 1) {
      s_ref <- sd(ref_units[, j])
      var_sum <- var_sum + s_ref^2 / n_ref
    }
    if (n_test > 1) {
      s_test <- sd(test_units[, j])
      var_sum <- var_sum + s_test^2 / n_test
    }
  }
  inner <- (1 / P) * (diff_sq + var_sum)
  if (!is.finite(inner) || inner < 0) return(NA_real_)
  100 - 25 * log10(1 + inner)
}

# Site/Python engine: 90% percentile CI of bootstrap E(f2); bc_f2 = 2*f2_obs - point E(f2).
python_style_bootstrap <- function(times, ref_units, test_units, n_iter = BOOT_ITER, seed = BOOT_SEED) {
  set.seed(seed)
  ref_mean <- colMeans(ref_units)
  test_mean <- colMeans(test_units)
  fda <- fda_truncate(times, ref_mean, test_mean)
  ema <- ema_truncate(times, ref_mean, test_mean)
  slices_diverge <- !identical(fda$idx, ema$idx)
  f2_obs <- if (slices_diverge) f2(ema$ref, ema$test) else f2(fda$ref, fda$test)
  point_e <- expected_f2_analytical(times, ref_units, test_units)
  n_ref <- nrow(ref_units)
  n_test <- nrow(test_units)
  boots <- numeric(n_iter)
  for (i in seq_len(n_iter)) {
    ri <- sample(seq_len(n_ref), n_ref, replace = TRUE)
    ti <- sample(seq_len(n_test), n_test, replace = TRUE)
    boots[i] <- expected_f2_analytical(
      times,
      ref_units[ri, , drop = FALSE],
      test_units[ti, , drop = FALSE]
    )
  }
  alpha <- 0.05
  ci <- quantile(boots, probs = c(alpha, 1 - alpha), type = 7)
  list(
    E_f2 = point_e,
    median = median(boots),
    bc_f2 = 2 * f2_obs - point_e,
    ci_lower = unname(ci[1]),
    ci_upper = unname(ci[2]),
    similar = unname(ci[1]) >= 50,
    f2_observed = f2_obs
  )
}

resolve_store_key <- function(ref_file, test_file, pair_name) {
  if (grepl("polli_slow_similar", ref_file, fixed = TRUE)) {
    return("polli_slow_similar")
  }
  if (grepl("mendyk", ref_file, fixed = TRUE)) {
    return("mendyk")
  }
  if (grepl("shah_ref", ref_file, fixed = TRUE)) {
    if (grepl("borderline", test_file, fixed = TRUE)) {
      return("shah_borderline")
    }
    return("shah")
  }
  if (grepl("tsong_ref", ref_file, fixed = TRUE)) {
    return("tsong")
  }
  if (grepl("chow_ref", ref_file, fixed = TRUE)) {
    return("chow")
  }
  pair_name
}

is_mendyk_pair <- function(ref_path, pair_label) {
  grepl("mendyk", basename(ref_path), fixed = TRUE) ||
    grepl("mendyk", pair_label, fixed = TRUE)
}

is_shah_pair <- function(ref_path, pair_label) {
  grepl("shah", basename(ref_path), fixed = TRUE) ||
    grepl("shah", pair_label, fixed = TRUE)
}

gold_scalar <- function(x) {
  if (is.null(x) || length(x) == 0) {
    return(NA_real_)
  }
  as.numeric(unlist(x))[1]
}

check_bootf2_shah1998_versions <- function() {
  cat("\n======== bootf2 version / shah1998 identity check ========\n")
  v_installed <- as.character(packageVersion("bootf2"))
  cat(sprintf("Installed bootf2 (active session): %s\n", v_installed))
  cat("Note: CRAN archive lists bootf2 0.4.0 only (no 1.0-2); comparing 0.4.1 vs archived 0.4.0.\n")

  ref_v041 <- as.data.frame(bootf2::shah1998$ref, stringsAsFactors = FALSE)
  test_v041 <- as.data.frame(bootf2::shah1998$test1, stringsAsFactors = FALSE)

  cmp_script <- tempfile(fileext = ".R")
  ref_out <- tempfile(fileext = ".rds")
  test_out <- tempfile(fileext = ".rds")
  on.exit({
    unlink(c(cmp_script, ref_out, test_out))
  }, add = TRUE)

  writeLines(
    c(
      sprintf('tmp_lib <- "%s"', gsub("\\\\", "/", file.path(tempdir(), "bootf2_cr040"))),
      "dir.create(tmp_lib, showWarnings=FALSE, recursive=TRUE)",
      'url <- "https://cran.r-project.org/src/contrib/Archive/bootf2/bootf2_0.4.0.tar.gz"',
      "install.packages(url, lib=tmp_lib, repos=NULL, type='source', quiet=TRUE)",
      "suppressPackageStartupMessages(library(bootf2, lib.loc=tmp_lib))",
      'cat("Archive bootf2:", as.character(packageVersion("bootf2")), "\\n")',
      "data(shah1998)",
      sprintf('saveRDS(shah1998$ref, "%s")', gsub("\\\\", "/", ref_out)),
      sprintf('saveRDS(shah1998$test1, "%s")', gsub("\\\\", "/", test_out))
    ),
    cmp_script
  )

  status <- system2("Rscript", cmp_script, stdout = TRUE, stderr = TRUE)
  if (length(status) > 0) {
    cat(paste(status, collapse = "\n"), "\n")
  }
  if (!file.exists(ref_out) || !file.exists(test_out)) {
    cat("*** WARNING: Could not load bootf2 0.4.0 for shah1998 comparison (install/build failed).\n")
    return(invisible(NULL))
  }

  ref_v040 <- readRDS(ref_out)
  test_v040 <- readRDS(test_out)
  ref_v040 <- as.data.frame(ref_v040, stringsAsFactors = FALSE)
  test_v040 <- as.data.frame(test_v040, stringsAsFactors = FALSE)

  ref_same <- isTRUE(all.equal(ref_v041, ref_v040, check.attributes = FALSE))
  test_same <- isTRUE(all.equal(test_v041, test_v040, check.attributes = FALSE))

  if (ref_same && test_same) {
    cat("shah1998$ref and shah1998$test1: BYTE-IDENTICAL across bootf2 0.4.1 and 0.4.0\n")
  } else {
    cat("*** ALERT: shah1998 data DIFFERS between bootf2 0.4.1 and 0.4.0 ***\n")
    if (!ref_same) {
      cat("  shah1998$ref: MISMATCH\n")
      print(all.equal(ref_v041, ref_v040, check.attributes = FALSE))
    }
    if (!test_same) {
      cat("  shah1998$test1: MISMATCH\n")
      print(all.equal(test_v041, test_v040, check.attributes = FALSE))
    }
  }
  invisible(list(ref_same = ref_same, test_same = test_same))
}

print_model_block <- function(summary, pair_label, side) {
  cat(sprintf(
    "\n--- %s | %s | %s | converged %d/%d ---\n",
    pair_label, side, summary$label, summary$n_converged, summary$n_total
  ))
  if (summary$n_converged == 0) {
    cat("  (no converged units)\n")
    return(invisible(NULL))
  }
  for (nm in names(summary$param_mean)) {
    cat(sprintf(
      "  %s = %.6f ± %.6f\n",
      nm, summary$param_mean[[nm]], summary$param_sd[[nm]]
    ))
  }
  g <- summary$gof_mean
  cat(sprintf(
    "  GOF: R2=%.4f, adj-R2=%.4f, RMSE=%.4f, AIC=%.4f, MSC=%.4f\n",
    g$R2, g$adj_R2, g$RMSE, g$AIC, g$MSC
  ))
}

run_pair_analysis <- function(ref_path, test_path, pair_label) {
  ref <- load_csv(ref_path)
  test <- load_csv(test_path)
  if (!identical(ref$times, test$times)) {
    stop(pair_label, ": reference and test time vectors must match")
  }

  cat(sprintf("\n######## %s ########\n", pair_label))

  ref_models <- list()
  for (mk in names(MODEL_REGISTRY)) {
    sm <- fit_model_all_units(ref$times, ref$units, mk)
    ref_models[[mk]] <- sm
    print_model_block(sm, pair_label, "ref")
  }

  metrics <- compute_profile_metrics(ref$times, ref$units, test$units)
  cat(sprintf("\n=== %s | similarity ===\n", pair_label))
  cat(sprintf(
    "  FDA:  f2=%.4f, f1=%.4f | times=%s\n",
    metrics$f2_fda, metrics$f1_fda,
    paste(format(metrics$truncated_times_fda, trim = TRUE), collapse = ",")
  ))
  cat(sprintf(
    "  EMA:  f2=%.4f, f1=%.4f | times=%s\n",
    metrics$f2_ema, metrics$f1_ema,
    paste(format(metrics$truncated_times_ema, trim = TRUE), collapse = ",")
  ))
  cat(sprintf(
    "  PRIMARY (%s): f2=%.4f, f1=%.4f\n",
    metrics$primary, metrics$f2_primary, metrics$f1_primary
  ))
  cat(sprintf(
    "  DE%% ref=%.4f, test=%.4f | MDT ref=%.4f, test=%.4f\n",
    metrics$DE_ref, metrics$DE_test, metrics$MDT_ref, metrics$MDT_test
  ))

  cv_ref <- cv_gate(ref$units)
  cv_test <- cv_gate(test$units)
  boot <- NULL
  if (is_shah_pair(ref_path, pair_label)) {
    cat(sprintf("\n=== %s | shah: deterministic f2 anchor — bootstrap skipped ===\n", pair_label))
  } else if (is_mendyk_pair(ref_path, pair_label) && (cv_ref$exceeded || cv_test$exceeded)) {
    cat(sprintf("\n=== %s | mendyk high-CV — bootstrap (EMA-primary / bootf2) ===\n", pair_label))
    if (cv_ref$exceeded) cat("  reference:", cv_ref$message, "\n")
    if (cv_test$exceeded) cat("  test:", cv_test$message, "\n")
    cat(sprintf(
      "  Seed: set.seed(%d) before bootf2::bootf2(..., n.boots=%d, regulation=\"EMA\")\n",
      BOOT_SEED, BOOT_ITER
    ))
    boot <- run_bootf2_ema(ref_path, test_path)
    py_boot <- python_style_bootstrap(ref$times, ref$units, test$units)
    cat(sprintf(
      "  bootf2 EMA: f2_obs=%.4f, est.f2 boot mean=%.4f, est.f2 boot median=%.4f, bc-f2=%.4f, bc-f2 90%% CI=[%.4f, %.4f], bc-f2 similar=%s\n",
      boot$f2_observed, boot$est_f2_boot_mean, boot$est_f2_boot_median, boot$bc_f2,
      boot$bc_f2_ci_lower, boot$bc_f2_ci_upper, boot$bc_f2_similar
    ))
    cat(sprintf(
      "  Python-aligned unit bootstrap (EMA, seed=%d): f2_obs=%.4f, E(f2)=%.4f, median=%.4f, bc-f2=%.4f, 90%% CI=[%.4f, %.4f], similar=%s\n",
      BOOT_SEED, py_boot$f2_observed, py_boot$E_f2, py_boot$median, py_boot$bc_f2,
      py_boot$ci_lower, py_boot$ci_upper, py_boot$similar
    ))
    boot$python_aligned <- py_boot
  } else if (is_mendyk_pair(ref_path, pair_label)) {
    cat(sprintf("\n=== %s | mendyk: CV gate OK — bootstrap not required ===\n", pair_label))
  } else {
    cat(sprintf("\n=== %s | CV gate OK — bootstrap not required ===\n", pair_label))
  }

  list(ref_models = ref_models, metrics = metrics, bootstrap = boot)
}

compare_row <- function(quantity, r_val, py_val, golden_val = NA, tol = NULL) {
  is_boot <- grepl("boot", quantity, ignore.case = TRUE)
  use_tol <- if (!is.null(tol)) tol else if (is_boot) 0.5 else 0.01
  r_val <- if (is.null(r_val) || length(r_val) == 0) NA_real_ else as.numeric(r_val)[1]
  py_val <- if (is.null(py_val) || length(py_val) == 0) NA_real_ else as.numeric(py_val)[1]
  golden_val <- if (length(golden_val) == 0 || is.null(golden_val)) NA_real_ else as.numeric(golden_val)[1]
  if (is.na(r_val) || is.na(py_val)) {
    status <- "MISSING"
    diff <- NA_real_
  } else {
    diff <- abs(r_val - py_val)
    status <- if (diff <= use_tol) "PASS" else "CHECK"
  }
  data.frame(
    quantity = quantity,
    R = ifelse(is.na(r_val), NA, round(r_val, 4)),
    Python = ifelse(is.na(py_val), NA, round(py_val, 4)),
    Golden = ifelse(is.na(golden_val), NA, round(golden_val, 4)),
    AbsDiff = ifelse(is.na(diff), NA, diff),
    Status = status,
    stringsAsFactors = FALSE
  )
}

run_validation <- function(results_by_pair, py_json_path, golden_path) {
  py_all <- jsonlite::fromJSON(py_json_path, simplifyVector = TRUE)
  golden <- if (!is.null(golden_path) && file.exists(golden_path)) {
    jsonlite::fromJSON(golden_path, simplifyVector = FALSE)
  } else {
    NULL
  }

  rows <- list()
  add <- function(...) {
    rows[[length(rows) + 1]] <<- compare_row(...)
  }

  add_ref_model <- function(prefix, r_model, py_block, golden_block = NULL) {
    if (is.null(r_model) || is.null(py_block)) {
      return(invisible(NULL))
    }
    for (p in names(r_model$param_mean)) {
      gval <- NA_real_
      if (!is.null(golden_block)) {
        gname <- if (p == "k1") "k" else p
        if (!is.null(golden_block[[gname]])) {
          gval <- gold_scalar(golden_block[[gname]])
        }
      }
      add(paste0(prefix, ".", p), r_model$param_mean[[p]], py_block[[p]], gval)
    }
    gof_map <- c(R2 = "R2", adj_R2 = "adjR2", RMSE = "RMSE", AIC = "AIC", MSC = "MSC")
    for (rk in names(gof_map)) {
      pyk <- gof_map[[rk]]
      gval <- NA_real_
      if (!is.null(golden_block) && !is.null(golden_block[[pyk]])) {
        gval <- gold_scalar(golden_block[[pyk]])
      } else if (!is.null(golden_block) && !is.null(golden_block[[rk]])) {
        gval <- gold_scalar(golden_block[[rk]])
      }
      add(
        paste0(prefix, ".", rk),
        r_model$gof_mean[[rk]],
        py_block[[pyk]],
        gval
      )
    }
  }

  add_similarity <- function(prefix, metrics, py, golden_sim = NULL) {
    add(paste0(prefix, ".f2_fda"), metrics$f2_fda, py$fda$f2,
        if (!is.null(golden_sim)) gold_scalar(golden_sim$f2) else NA)
    add(paste0(prefix, ".f1_fda"), metrics$f1_fda, py$fda$f1,
        if (!is.null(golden_sim)) gold_scalar(golden_sim$f1) else NA)
    add(paste0(prefix, ".f2_ema"), metrics$f2_ema, py$ema$f2, NA)
    add(paste0(prefix, ".f1_ema"), metrics$f1_ema, py$ema$f1, NA)
    add(paste0(prefix, ".DE_ref"), metrics$DE_ref, py$fda$DE_ref,
        if (!is.null(golden_sim)) gold_scalar(golden_sim$DE_ref) else NA)
    add(paste0(prefix, ".DE_test"), metrics$DE_test, py$fda$DE_test,
        if (!is.null(golden_sim)) gold_scalar(golden_sim$DE_test) else NA)
    add(paste0(prefix, ".MDT_ref"), metrics$MDT_ref, py$fda$MDT_ref,
        if (!is.null(golden_sim)) gold_scalar(golden_sim$MDT_ref) else NA)
    add(paste0(prefix, ".MDT_test"), metrics$MDT_test, py$fda$MDT_test,
        if (!is.null(golden_sim)) gold_scalar(golden_sim$MDT_test) else NA)
  }

  if (!is.null(results_by_pair$polli_slow_similar)) {
    r <- results_by_pair$polli_slow_similar
    py <- py_all$polli_slow
    gsec <- golden$polli_slow_similar
    gfits <- gsec$reference_model_fits
    add_ref_model("polli_slow.Weibull", r$ref_models$weibull, py$weibull, gfits$Weibull)
    add_ref_model("polli_slow.FO", r$ref_models$first_order, py$fo, gfits[["First-order"]])
    add_ref_model("polli_slow.HC", r$ref_models$hixson_crowell, py$hc, gfits[["Hixson-Crowell"]])
    add_ref_model("polli_slow.KP", r$ref_models$korsmeyer_peppas_full, py$kp,
                  gfits[["Korsmeyer-Peppas_full"]])
    add_similarity("polli_slow", r$metrics, py, gsec$similarity_ref_vs_test_unit6)
  }

  if (!is.null(results_by_pair$mendyk)) {
    r <- results_by_pair$mendyk
    py <- py_all$mendyk
    add_ref_model("mendyk.Weibull", r$ref_models$weibull, py$weibull)
    add_ref_model("mendyk.FO", r$ref_models$first_order, py$fo)
    add_ref_model("mendyk.HC", r$ref_models$hixson_crowell, py$hc)
    add_ref_model("mendyk.KP", r$ref_models$korsmeyer_peppas_full, py$kp)
    add_similarity("mendyk", r$metrics, py)
    if (!is.null(r$bootstrap$python_aligned)) {
      py_boot <- py$boot_ema
      pb <- r$bootstrap$python_aligned
      bf <- r$bootstrap
      add("mendyk.boot.f2_observed", pb$f2_observed, py_boot$f2_observed, NA, 0.01)
      add("mendyk.boot.E_f2", pb$E_f2, py_boot$E_f2, NA, 0.5)
      add("mendyk.boot.median", pb$median, py_boot$median, NA, 0.5)
      add("mendyk.boot.bc_f2", pb$bc_f2, py_boot$bc_f2, NA, 0.5)
      add("mendyk.boot.ci_lower", pb$ci_lower, py_boot$ci_lower, NA, 1.5)
      add("mendyk.boot.ci_upper", pb$ci_upper, py_boot$ci_upper, NA, 1.5)
    }
  }

  if (!is.null(results_by_pair$shah)) {
    r <- results_by_pair$shah
    py <- py_all$shah
    add_ref_model("shah.Weibull", r$ref_models$weibull, py$weibull)
    add_ref_model("shah.FO", r$ref_models$first_order, py$fo)
    add_ref_model("shah.HC", r$ref_models$hixson_crowell, py$hc)
    add_ref_model("shah.KP", r$ref_models$korsmeyer_peppas_full, py$kp)
    add_similarity("shah", r$metrics, py)
  }

  if (!is.null(results_by_pair$shah_borderline)) {
    r <- results_by_pair$shah_borderline
    py <- py_all$shah_borderline
    add_ref_model("shah_borderline.Weibull", r$ref_models$weibull, py$weibull)
    add_ref_model("shah_borderline.FO", r$ref_models$first_order, py$fo)
    add_ref_model("shah_borderline.HC", r$ref_models$hixson_crowell, py$hc)
    add_ref_model("shah_borderline.KP", r$ref_models$korsmeyer_peppas_full, py$kp)
    add_similarity("shah_borderline", r$metrics, py)
  }

  if (!is.null(results_by_pair$tsong)) {
    r <- results_by_pair$tsong
    py <- py_all$tsong
    add_ref_model("tsong.Weibull", r$ref_models$weibull, py$weibull)
    add_ref_model("tsong.FO", r$ref_models$first_order, py$fo)
    add_ref_model("tsong.HC", r$ref_models$hixson_crowell, py$hc)
    add_ref_model("tsong.KP", r$ref_models$korsmeyer_peppas_full, py$kp)
    add_similarity("tsong", r$metrics, py)
  }

  if (!is.null(results_by_pair$chow)) {
    r <- results_by_pair$chow
    py <- py_all$chow
    add_ref_model("chow.Weibull", r$ref_models$weibull, py$weibull)
    add_ref_model("chow.FO", r$ref_models$first_order, py$fo)
    add_ref_model("chow.HC", r$ref_models$hixson_crowell, py$hc)
    add_ref_model("chow.KP", r$ref_models$korsmeyer_peppas_full, py$kp)
    add_similarity("chow", r$metrics, py)
  }

  tab <- do.call(rbind, rows)
  n_pass <- sum(tab$Status == "PASS", na.rm = TRUE)
  n_check <- sum(tab$Status == "CHECK", na.rm = TRUE)
  n_miss <- sum(tab$Status == "MISSING", na.rm = TRUE)
  cat("\n\n======== VALIDATION: R vs Python (+ golden for polli_slow_similar) ========\n")
  cat(sprintf(
    "Packages: minpack.lm %s | bootf2 %s | jsonlite %s\n",
    as.character(packageVersion("minpack.lm")),
    as.character(packageVersion("bootf2")),
    as.character(packageVersion("jsonlite"))
  ))
  cat(sprintf("Bootstrap MC tolerance: E(f2)/bc-f2 ±0.5; CI endpoints ±1.5; deterministic ±0.01\n"))
  print(tab, row.names = FALSE)
  cat(sprintf(
    "\nValidation summary: %d PASS, %d CHECK, %d MISSING\n",
    n_pass, n_check, n_miss
  ))
  boot_rows <- tab[grepl("^mendyk\\.boot\\.", tab$quantity) & !grepl("bootf2", tab$quantity), ]
  if (nrow(boot_rows) > 0 && all(boot_rows$Status == "PASS", na.rm = TRUE)) {
    cat("Mendyk bootstrap (Python-aligned EMA-primary): agreement within Monte-Carlo tolerance.\n")
  }
  invisible(tab)
}

main <- function() {
  if (!requireNamespace("minpack.lm", quietly = TRUE)) {
    stop("Install minpack.lm: install.packages('minpack.lm')")
  }
  if (!requireNamespace("bootf2", quietly = TRUE)) {
    stop("Install bootf2: install.packages('bootf2')")
  }
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    stop("Install jsonlite: install.packages('jsonlite')")
  }

  check_bootf2_shah1998_versions()

  cfg <- parse_args(args)
  pair_specs <- lapply(cfg$pairs, function(p) {
    parts <- strsplit(trimws(p), ":", fixed = TRUE)[[1]]
    list(ref = parts[1], test = parts[2])
  })

  results_by_pair <- list()
  for (spec in pair_specs) {
    ref_path <- file.path(cfg$data_dir, spec$ref)
    test_path <- file.path(cfg$data_dir, spec$test)
    pair_name <- paste0(
      tools::file_path_sans_ext(spec$ref), "_vs_",
      tools::file_path_sans_ext(spec$test)
    )
    if (grepl("polli_slow_similar", spec$ref)) {
      store_key <- "polli_slow_similar"
    } else {
      store_key <- resolve_store_key(spec$ref, spec$test, pair_name)
    }
    results_by_pair[[store_key]] <- run_pair_analysis(ref_path, test_path, pair_name)
  }

  if (!is.null(cfg$validate_json)) {
    run_validation(results_by_pair, cfg$validate_json, cfg$golden_json)
  }
}

main()
