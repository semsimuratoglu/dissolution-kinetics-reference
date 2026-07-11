# Canonical Shah datasets exported from bootf2::shah1998:
#   shah_ref.csv              <- shah1998$ref
#   shah_test.csv             <- shah1998$test1   (shah_similar, f2≈60.03)
#   shah_borderline_test.csv  <- shah1998$test5   (shah_borderline, f2≈48.05)
#
# Usage:
#   Rscript scripts/export_shah1998.R [--out-dir data]

args <- commandArgs(trailingOnly = TRUE)

parse_args <- function(args) {
  out_dir <- "data"
  i <- 1
  while (i <= length(args)) {
    if (args[i] == "--out-dir" && i < length(args)) {
      out_dir <- args[i + 1]
      i <- i + 2
    } else {
      i <- i + 1
    }
  }
  list(out_dir = out_dir)
}

bootf2_df_to_csv <- function(df, path) {
  if (!"time" %in% names(df)) {
    stop("Expected a 'time' column in bootf2 data frame")
  }
  unit_cols <- setdiff(names(df), "time")
  n_units <- length(unit_cols)
  out <- data.frame(Time = df$time)
  for (j in seq_len(n_units)) {
    out[[paste0("Rep ", j)]] <- df[[unit_cols[j]]]
  }
  write.csv(out, path, row.names = FALSE, quote = FALSE)
  invisible(path)
}

cfg <- parse_args(args)

if (!requireNamespace("bootf2", quietly = TRUE)) {
  stop("Install bootf2: install.packages('bootf2')")
}
if (!requireNamespace("jsonlite", quietly = TRUE)) {
  stop("Install jsonlite: install.packages('jsonlite')")
}

bootf2_version <- as.character(packageVersion("bootf2"))
data(shah1998, package = "bootf2", envir = environment())

out_dir <- cfg$out_dir
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

ref_path <- file.path(out_dir, "shah_ref.csv")
test_path <- file.path(out_dir, "shah_test.csv")
borderline_path <- file.path(out_dir, "shah_borderline_test.csv")

bootf2_df_to_csv(shah1998$ref, ref_path)
bootf2_df_to_csv(shah1998$test1, test_path)
bootf2_df_to_csv(shah1998$test5, borderline_path)

cat(sprintf("bootf2 version: %s\n", bootf2_version))
cat(sprintf("Exported shah_ref.csv              <- shah1998$ref   (%d time points, %d units)\n",
            nrow(shah1998$ref), ncol(shah1998$ref) - 1))
cat(sprintf("Exported shah_test.csv             <- shah1998$test1  (%d time points, %d units)\n",
            nrow(shah1998$test1), ncol(shah1998$test1) - 1))
cat(sprintf("Exported shah_borderline_test.csv  <- shah1998$test5  (%d time points, %d units)\n",
            nrow(shah1998$test5), ncol(shah1998$test5) - 1))
cat(sprintf("  -> %s\n", ref_path))
cat(sprintf("  -> %s\n", test_path))
cat(sprintf("  -> %s\n", borderline_path))

meta <- list(
  source = "bootf2::shah1998",
  bootf2_version = bootf2_version,
  citation = "Shah et al., Pharm. Res. 1998, 15(6), 889-896",
  pairs = list(
    shah_similar = list(
      reference_arm = "ref",
      test_arm = "test1",
      ref_file = "shah_ref.csv",
      test_file = "shah_test.csv",
      role = "deterministic f2 anchor (≈60.03); bootstrap not required"
    ),
    shah_borderline = list(
      reference_arm = "ref",
      test_arm = "test5",
      ref_file = "shah_ref.csv",
      test_file = "shah_borderline_test.csv",
      role = "borderline f2 (≈48.05); bootstrap not required under engine CV gate"
    )
  ),
  time_points_min = as.numeric(shah1998$ref$time),
  n_units = ncol(shah1998$ref) - 1
)
meta_path <- file.path(out_dir, "shah1998_provenance.json")
writeLines(jsonlite::toJSON(meta, pretty = TRUE, auto_unbox = TRUE), meta_path)
cat(sprintf("Provenance: %s\n", meta_path))
