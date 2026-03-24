library(LDZipMatrix)
library(tinytest)

source("../../../scripts/source.R")

# ============================================================================
# Test 1: Basic pairwise functionality with integer indices
# ============================================================================

cat("\n=== Testing basic pairwise with integer indices ===\n")

# Load test data - use lossless compression
orig <- read_plink_ld_matrix("../../tests/data/test_index.bin")
ld <- LDZipMatrix("../../tests/data/compressed_bits_index_99_min_0")

# Test single pair
result <- fetchLD(ld, 1, 1, types="PHASED_R", pairwise=TRUE)
expect_true(is.data.frame(result), info="Should return data.frame")
expect_equal(nrow(result), 1, info="Should have 1 row")
expect_equal(result$PHASED_R[1], orig[1,1], info="Value should match")

# Test multiple pairs - diagonal
result <- fetchLD(ld, 1:5, 1:5, types="PHASED_R", pairwise=TRUE)
expect_equal(nrow(result), 5, info="Should have 5 rows")
for (i in 1:5) {
  expect_equal(result$PHASED_R[i], orig[i,i], tolerance=1e-7,
               info=sprintf("Diagonal value (%d,%d) should match", i, i))
}

# Test multiple pairs - off-diagonal
rows <- c(1, 5, 10)
cols <- c(2, 8, 15)
result <- fetchLD(ld, rows, cols, types="PHASED_R", pairwise=TRUE)
expect_equal(nrow(result), 3, info="Should have 3 rows")
# Results are in input order
for (i in seq_along(rows)) {
  expect_equal(result$PHASED_R[i], orig[rows[i], cols[i]], tolerance=1e-7,
               info=sprintf("Value (%d,%d) should match", rows[i], cols[i]))
  expect_equal(result$var1[i], rows[i], info="var1 should match input row")
  expect_equal(result$var2[i], cols[i], info="var2 should match input col")
}

# Test same column, different rows
rows <- c(1, 5, 10)
cols <- c(50, 50, 50)
result <- fetchLD(ld, rows, cols, types="PHASED_R", pairwise=TRUE)
expect_equal(nrow(result), 3, info="Should have 3 rows")
# Results are in input order
for (i in seq_along(rows)) {
  expect_equal(result$PHASED_R[i], orig[rows[i], cols[i]], tolerance=1e-7,
               info=sprintf("Value (%d,%d) should match", rows[i], cols[i]))
}

# ============================================================================
# Test 2: Pairwise with rsids
# ============================================================================

cat("\n=== Testing pairwise with rsids ===\n")

orig_rsid <- read_tabular_ld_matrix("../../tests/data/test_tabular")
ld_rsid <- LDZipMatrix("../../tests/data/compressed_bits_tabular_99_min_0")
suppressMessages(buildIndex(ld_rsid))

# Test with rsids
test_rsids_row <- c("rs1", "rs10", "rs50")
test_rsids_col <- c("rs5", "rs15", "rs55")

result <- fetchLD(ld_rsid, test_rsids_row, test_rsids_col, types="UNPHASED_R", pairwise=TRUE)
expect_true(is.data.frame(result), info="Result should be data.frame")
expect_equal(nrow(result), 3, info="Result should have 3 rows")

# Verify values (results are in input order)
for (i in seq_along(test_rsids_row)) {
  expect_equal(result$UNPHASED_R[i], orig_rsid[test_rsids_row[i], test_rsids_col[i]], tolerance=1e-7,
               info=sprintf("Value (%s,%s) should match", test_rsids_row[i], test_rsids_col[i]))
  expect_equal(result$var1[i], test_rsids_row[i], info="var1 should match input")
  expect_equal(result$var2[i], test_rsids_col[i], info="var2 should match input")
}

# ============================================================================
# Test 3: Data frame structure tests
# ============================================================================

cat("\n=== Testing data frame structure ===\n")

# Single statistic → data.frame with stat + var1 + var2 columns
result <- fetchLD(ld, c(1,5,10), c(2,8,15), types="PHASED_R", pairwise=TRUE)
expect_true(is.data.frame(result), info="Pairwise should return data.frame")
expect_equal(nrow(result), 3, info="Result should have 3 rows")
expect_equal(ncol(result), 3, info="Single stat should give 3 columns (stat + var1 + var2)")
expect_true("PHASED_R" %in% colnames(result), info="Column name should be PHASED_R")
expect_true("var1" %in% colnames(result), info="Should have var1 column")
expect_true("var2" %in% colnames(result), info="Should have var2 column")

# ============================================================================
# Test 4: Error conditions
# ============================================================================

cat("\n=== Testing error conditions ===\n")

# Test mismatched lengths
expect_error(
  fetchLD(ld, c(1, 2, 3), c(1, 2), types="PHASED_R", pairwise=TRUE),
  pattern = "equal length",
  info = "Should fail when cols and rows have different lengths"
)

# Test invalid stat type
expect_error(
  fetchLD(ld, c(1, 2), c(1, 2), types="INVALID_STAT", pairwise=TRUE),
  pattern = "Unsupported Stat",
  info = "Should fail with invalid stat type"
)

# Test that for a given column, rows must be sorted (C++ requirement)
# Pairs: (col=5, row=10), (col=5, row=2), (col=5, row=8)
# All have same column, but rows are unsorted -> should error from C++
expect_error(
  fetchLD(ld, c(10, 2, 8), c(5, 5, 5), types="PHASED_R", pairwise=TRUE),
  pattern = "rows.*sorted",
  info = "Should fail when rows for same column are not sorted"
)

# Test that for different columns with unsorted rows works if each column's rows are sorted
# Pairs: (col=5, row=1), (col=5, row=2), (col=10, row=3), (col=10, row=4)
# col=5 has rows [1,2] (sorted), col=10 has rows [3,4] (sorted) -> should work
result_sorted_per_col <- fetchLD(ld, c(1, 2, 3, 4), c(5, 5, 10, 10), types="PHASED_R", pairwise=TRUE)
expect_equal(nrow(result_sorted_per_col), 4, info="Should work when rows are sorted per column")

# Test duplicate pairs are now allowed (inputs must be sorted)
result_dup <- fetchLD(ld, c(1, 1, 2), c(2, 2, 3), types="PHASED_R", pairwise=TRUE)
expect_true(is.data.frame(result_dup), info="Should succeed with duplicate pairs")
expect_equal(nrow(result_dup), 3, info="Should return 3 rows including duplicates")
expect_equal(result_dup$PHASED_R[1], result_dup$PHASED_R[2], info="Duplicate pairs should have same value")
expect_equal(result_dup$PHASED_R[1], orig[1,2], tolerance=1e-7, info="First duplicate should match original")
expect_equal(result_dup$PHASED_R[2], orig[1,2], tolerance=1e-7, info="Second duplicate should match original")
expect_equal(result_dup$var1[1], 1, info="First duplicate var1 should be 1")
expect_equal(result_dup$var2[1], 2, info="First duplicate var2 should be 2")
expect_equal(result_dup$var1[2], 1, info="Second duplicate var1 should be 1")
expect_equal(result_dup$var2[2], 2, info="Second duplicate var2 should be 2")

result_dup2 <- fetchLD(ld, c(5, 5, 10), c(10, 10, 20), types="PHASED_R", pairwise=TRUE)
expect_true(is.data.frame(result_dup2), info="Should succeed with duplicate pairs")
expect_equal(nrow(result_dup2), 3, info="Should return 3 rows including duplicates")
expect_equal(result_dup2$PHASED_R[1], result_dup2$PHASED_R[2], info="Duplicate pairs should have same value")

