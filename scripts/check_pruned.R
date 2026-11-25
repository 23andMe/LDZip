#!/usr/bin/env Rscript

# Check LD pruning results by validating against original vcor file
# Usage: Rscript check_pruned.R <vcor_file> <vars_file> <prune_in> <prune_out> <threshold> <window_kb>

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 6) {
  stop("Usage: Rscript check_pruned.R <vcor_file> <vars_file> <prune_in> <prune_out> <threshold> <window_kb>")
}

vcor_file <- args[1]
vars_file <- args[2]
prune_in_file <- args[3]
prune_out_file <- args[4]
threshold <- as.numeric(args[5])
window_kb <- as.numeric(args[6])

cat("Reading files...\n")

# Read LD data
ld_data <- read.table(vcor_file, header = TRUE, sep = "\t", comment.char = "", check.names = FALSE)
cat(sprintf("  LD pairs: %d\n", nrow(ld_data)))

# Read variant positions
vars <- read.table(vars_file, header = TRUE, sep = "\t", comment.char = "")
cat(sprintf("  Variants: %d\n", nrow(vars)))

# Read pruning results (0-based indices)
kept <- as.integer(readLines(prune_in_file))
removed <- as.integer(readLines(prune_out_file))
cat(sprintf("  Kept: %d, Removed: %d\n", length(kept), length(removed)))

# Get LD values from last column
ld_col <- ncol(ld_data)
ld_data$r_val <- ld_data[[ld_col]]

# Map variant IDs to indices
id_to_idx <- setNames(0:(nrow(vars) - 1), vars$ID)
ld_data$idx_A <- id_to_idx[as.character(ld_data$`#ID_A`)]
ld_data$idx_B <- id_to_idx[as.character(ld_data$ID_B)]
ld_data <- ld_data[!is.na(ld_data$idx_A) & !is.na(ld_data$idx_B), ]

# Add distance column
ld_data$dist_bp <- abs(vars$POS[ld_data$idx_A + 1] - vars$POS[ld_data$idx_B + 1])

cat("\nValidating pruning results...\n")

# Test: For high LD pairs within window, both should NOT be in kept
cat("\n[Test] Checking high-LD pairs from vcor file...\n")

# Filter to pairs within window
in_window <- ld_data$dist_bp <= window_kb * 1000
window_pairs <- ld_data[in_window, ]
cat(sprintf("  Pairs within window: %d\n", nrow(window_pairs)))

# Filter to high LD
high_ld <- abs(window_pairs$r_val) >= threshold
high_ld_pairs <- window_pairs[high_ld, ]
cat(sprintf("  High-LD pairs (|r| >= %.4f): %d\n", threshold, nrow(high_ld_pairs)))

if (nrow(high_ld_pairs) > 0) {
  # Sample random pairs to check
  n_tests <- min(1000, nrow(high_ld_pairs))
  test_rows <- sample(1:nrow(high_ld_pairs), n_tests)

  violations <- 0
  for (row_idx in test_rows) {
    row <- high_ld_pairs[row_idx, ]
    idx_A <- row$idx_A
    idx_B <- row$idx_B

    both_kept <- (idx_A %in% kept) && (idx_B %in% kept)

    if (both_kept) {
      violations <- violations + 1
      if (violations <= 5) {
        cat(sprintf("  ❌ Violation: variants %d and %d both kept with |r|=%.4f (dist=%d bp)\n",
                    idx_A, idx_B, abs(row$r_val), row$dist_bp))
      }
    }
  }

  cat(sprintf("  Tested %d random high-LD pairs\n", n_tests))
  if (violations > 0) {
    cat(sprintf("  ❌ Found %d violations: both variants kept despite high LD\n", violations))
    quit(status = 1)
  } else {
    cat("  ✅ All high-LD pairs: at least one variant removed\n")
  }
} else {
  cat("  ℹ️  No high-LD pairs within window found\n")
}

cat("\n✅ Validation complete!\n")
