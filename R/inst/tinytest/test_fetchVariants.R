library(LDZipMatrix)
library(tinytest)

# Test with tabular data (has rsIDs and indexed)
bits_list <- c(8, 16, 32, 99)
min_list  <- c(0.00, 0.001, 0.01, 0.1)

for (bits in bits_list) {
  for (min in min_list) {
    prefix <- sprintf("../../tests/data/compressed_bits_tabular_%d_min_%s", bits, min)
    cat("Testing:", prefix, "\n")

    # load compressed
    ld <- LDZipMatrix(prefix)
    suppressMessages(buildIndex(ld))

    # Test 1: Query by numeric indices
    indices <- c(1, 10, 50, 100)
    result_idx <- fetchVariants(ld, indices)

    expect_inherits(result_idx, "data.frame",
                    info = sprintf("fetchVariants should return data.frame for %s", prefix))

    expect_equal(nrow(result_idx), length(indices),
                 info = sprintf("Should return %d rows for %d indices", length(indices), length(indices)))

    expect_true(all(c("idx", "CHROM", "POS", "ID", "REF", "ALT") %in% colnames(result_idx)),
                info = sprintf("Result should have all required columns for %s", prefix))

    expect_equal(result_idx$idx, indices,
                 info = sprintf("Indices should match input for %s", prefix))

    # Test 2: Query by rsIDs
    rsids <- c("rs1", "rs10", "rs50", "rs100")
    result_rsid <- fetchVariants(ld, rsids)

    expect_equal(nrow(result_rsid), length(rsids),
                 info = sprintf("Should return %d rows for %d rsIDs", length(rsids), length(rsids)))

    expect_equal(result_rsid$ID, rsids,
                 info = sprintf("IDs should match input rsIDs for %s", prefix))

    expect_equal(result_rsid$idx, indices,
                 info = sprintf("rsID query should return same indices as numeric query for %s", prefix))

    # Test 3: Query by genomic region
    region <- "2:10000-10020"
    result_region <- fetchVariants(ld, region)

    expect_true(nrow(result_region) > 0,
                info = sprintf("Region query should return results for %s", prefix))

    expect_equal(result_region$CHROM, rep(2, nrow(result_region)),
                 info = sprintf("All variants should be on chromosome 2 for %s", prefix))

    expect_true(all(result_region$POS >= 10000 & result_region$POS <= 10020),
                info = sprintf("All positions should be within region bounds for %s", prefix))

    # Expected variants: rs1 to rs21 (POS 10000 to 10020, inclusive)
    expected_rsids <- paste0("rs", 1:21)
    expect_equal(result_region$ID, expected_rsids,
                 info = sprintf("Region should return rs1-rs21 for %s", prefix))

    # Test 4: Query with "chr" prefix in region
    region_chr <- "chr2:10000-10020"
    result_region_chr <- fetchVariants(ld, region_chr)

    expect_equal(result_region_chr, result_region,
                info = sprintf("Region with 'chr' prefix should match without for %s", prefix))

    # Test 5: Single variant queries
    single_idx <- fetchVariants(ld, 5)
    expect_equal(nrow(single_idx), 1,
                 info = sprintf("Single index should return 1 row for %s", prefix))
    expect_equal(single_idx$idx, 5,
                 info = sprintf("Single index query should return correct idx for %s", prefix))

    single_rsid <- fetchVariants(ld, "rs5")
    expect_equal(nrow(single_rsid), 1,
                 info = sprintf("Single rsID should return 1 row for %s", prefix))
    expect_equal(single_rsid$ID, "rs5",
                 info = sprintf("Single rsID query should return correct ID for %s", prefix))

    # Test 6: Results are sorted by index
    unsorted_indices <- c(100, 10, 50, 1, 25)
    result_unsorted <- fetchVariants(ld, unsorted_indices)
    expect_true(all(diff(result_unsorted$idx) > 0),
                info = sprintf("Results should be sorted by idx for %s", prefix))
    expect_equal(result_unsorted$idx, sort(unsorted_indices),
                 info = sprintf("Sorted indices should match sorted input for %s", prefix))

    # Test 7: Verify CPRA completeness (no NAs in critical fields)
    result_all <- fetchVariants(ld, 1:20)
    expect_true(all(!is.na(result_all$CHROM)),
                info = sprintf("CHROM should not have NAs for %s", prefix))
    expect_true(all(!is.na(result_all$POS)),
                info = sprintf("POS should not have NAs for %s", prefix))
    expect_true(all(!is.na(result_all$ID)),
                info = sprintf("ID should not have NAs for %s", prefix))

    # Test 8: Verify REF and ALT are present (may be NA depending on data)
    expect_true("REF" %in% colnames(result_all),
                info = sprintf("REF column should exist for %s", prefix))
    expect_true("ALT" %in% colnames(result_all),
                info = sprintf("ALT column should exist for %s", prefix))
  }
}

# Test error conditions with a single prefix
prefix <- "../../tests/data/compressed_bits_tabular_16_min_0"
ld <- LDZipMatrix(prefix)
suppressMessages(buildIndex(ld))

# Test 9: Empty input should error
expect_error(
  fetchVariants(ld, integer(0)),
  "`variants` must be non-empty",
  info = "Empty input should error"
)

expect_error(
  fetchVariants(ld, character(0)),
  "`variants` must be non-empty",
  info = "Empty character vector should error"
)

# Test 10: Duplicate entries should error
expect_error(
  fetchVariants(ld, c(1, 1)),
  "duplicate",
  info = "Duplicate indices should error"
)

expect_error(
  fetchVariants(ld, c("rs1", "rs1")),
  "duplicate",
  info = "Duplicate rsIDs should error"
)

# Test 11: Invalid numeric indices
expect_error(
  fetchVariants(ld, 0),
  "positive",
  info = "Index zero should error"
)

expect_error(
  fetchVariants(ld, 1.5),
  "integer",
  info = "Non-integer index should error"
)

# Test 13: Out-of-range region should error
expect_error(
  fetchVariants(ld, "2:99000000-99999999"),
  "No matching variants found in specified region",
  info = "Out-of-range region should error"
)

# Test 14: Non-existent or invalid rsID should error
expect_error(
  fetchVariants(ld, c("rs1", "rs_nonexistent")),
  "No numeric part found in identifier",
  info = "Invalid rsID format should error"
)

expect_error(
  fetchVariants(ld, "rs99999999"),
  "No matching rsIDs found",
  info = "Non-existent rsID should error"
)

# Test 15: Missing database file (no buildIndex)
prefix_no_index <- "../../tests/data/compressed_bits_index_16_min_0"
ld_no_index <- LDZipMatrix(prefix_no_index)

expect_error(
  fetchVariants(ld_no_index, "rs1"),
  "Database file not found.*buildIndex",
  info = "Missing database should give helpful error"
)
