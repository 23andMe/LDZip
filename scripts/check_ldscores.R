#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) {
  stop("Usage: Rscript check_ldscores.R <vcor_file> <vars_file> <ldscores_file> <window_kb>")
}

vcor_file <- args[1]
vars_file <- args[2]
ldscores_file <- args[3]
window_kb <- as.numeric(args[4])

cat("Reading files...\n")

ld_data <- read.table(vcor_file, header = TRUE, sep = "\t", comment.char = "", check.names = FALSE)
cat(sprintf("  LD pairs: %d\n", nrow(ld_data)))

vars <- read.table(vars_file, header = TRUE, sep = "\t", comment.char = "")
cat(sprintf("  Variants: %d\n", nrow(vars)))

ldscores <- read.table(ldscores_file, header = TRUE, sep = "\t")
cat(sprintf("  LD scores: %d\n", nrow(ldscores)))

if (nrow(ldscores) != nrow(vars)) {
  stop("LD score count does not match variant count")
}

id_to_idx <- setNames(0:(nrow(vars) - 1), vars$ID)
ld_data$idx_A <- id_to_idx[as.character(ld_data$`#ID_A`)]
ld_data$idx_B <- id_to_idx[as.character(ld_data$ID_B)]
ld_data <- ld_data[!is.na(ld_data$idx_A) & !is.na(ld_data$idx_B), ]

ld_col <- ncol(ld_data) - 2
ld_data$r2_val <- ld_data[[ld_col]]

ld_data$dist_bp <- abs(vars$POS[ld_data$idx_A + 1] - vars$POS[ld_data$idx_B + 1])

cat("\nCalculating expected LD scores...\n")

window_bp <- window_kb * 1000
n_variants <- nrow(vars)
expected_ldscores <- numeric(n_variants)

for (i in 1:nrow(ld_data)) {
  if (ld_data$dist_bp[i] > window_bp) next

  idx_A <- ld_data$idx_A[i]
  idx_B <- ld_data$idx_B[i]
  r2 <- ld_data$r2_val[i]

  expected_ldscores[idx_A + 1] <- expected_ldscores[idx_A + 1] + r2
  expected_ldscores[idx_B + 1] <- expected_ldscores[idx_B + 1] + r2
}

cat("\nComparing LD scores...\n")

observed <- ldscores$ld_score
expected <- expected_ldscores

n_test <- min(100, n_variants)
test_indices <- sample(1:n_variants, n_test)

max_diff <- 0
violations <- 0

tolerance <- 5e-4

for (i in test_indices) {
  diff <- abs(observed[i] - expected[i])
  max_diff <- max(max_diff, diff)

  if (diff > tolerance) {
    violations <- violations + 1
    if (violations <= 5) {
      cat(sprintf("  Variant %d: observed=%.6f, expected=%.6f, diff=%.6f\n",
                  i - 1, observed[i], expected[i], diff))
    }
  }
}

cat(sprintf("\nTested %d variants\n", n_test))
cat(sprintf("Max difference: %.10f\n", max_diff))

if (violations > 0) {
  cat(sprintf("❌ Found %d variants with differences > %g\n", violations, tolerance))
  quit(status = 1)
} else {
  cat("✅ All LD scores match within tolerance!\n")
}
