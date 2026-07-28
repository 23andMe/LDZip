library(LDZipMatrix)
library(tinytest)

# Test all examples from README.md
prefix <- file.path(system.file("extdata", package = "LDZipMatrix"), "g1k.chr22.ldzip")

# --- LDZipMatrix() constructor ---
ld <- LDZipMatrix(prefix)
expect_inherits(ld, "LDZipMatrix", info = "LDZipMatrix() should return LDZipMatrix object")

# --- buildIndex() ---
buildIndex(ld)
db_file <- paste0(prefix, ".sqlite")
expect_true(file.exists(db_file), info = "buildIndex() should create .sqlite file")

# --- fetchLD() examples ---

# Example 1: Query by integer indices with multiple types
result1 <- fetchLD(ld, 2, 191, types = c("PHASED_R", "DPRIME"))
expect_true(!is.null(result1), info = "fetchLD with indices and multiple types should work")
if (is.data.frame(result1)) {
  expect_true(all(c("PHASED_R", "DPRIME") %in% colnames(result1)),
              info = "Should have both PHASED_R and DPRIME columns")
} else if (is.list(result1)) {
  expect_equal(names(result1), c("PHASED_R", "DPRIME"),
               info = "List should have PHASED_R and DPRIME elements")
}

# Example 2: Query by rsIDs
result2 <- fetchLD(ld, "rs587755077", "rs587631919", types = "PHASED_R")
expect_true(is.numeric(result2) || is.data.frame(result2) || is.list(result2),
            info = "fetchLD with rsIDs should return numeric/data.frame/list")

# Example 3: Query by region
result3 <- fetchLD(ld, "22:16050000-16051000", "22:16050000-16051000", types = "PHASED_R")
expect_true(is.matrix(result3) || is.list(result3),
            info = "fetchLD with region should return matrix or list")

# --- fetchVariants() examples ---

# Example 1: Query by numeric indices
v1 <- fetchVariants(ld, c(1, 10, 50))
expect_inherits(v1, "data.frame", info = "fetchVariants with indices should return data.frame")
expect_equal(nrow(v1), 3, info = "Should return 3 rows")
expect_equal(v1$idx, c(1, 10, 50), info = "Indices should match input")
expect_true(all(c("idx", "CHROM", "POS", "ID", "REF", "ALT") %in% colnames(v1)),
            info = "Should have all CPRA columns")

# Example 2: Query by rsIDs
v2 <- fetchVariants(ld, c("rs587755077", "rs587631919"))
expect_inherits(v2, "data.frame", info = "fetchVariants with rsIDs should return data.frame")
expect_equal(nrow(v2), 2, info = "Should return 2 rows")
expect_equal(v2$ID, c("rs587755077", "rs587631919"), info = "IDs should match input")

# Example 3: Query by region
v3 <- fetchVariants(ld, "22:16050000-16051000")
expect_inherits(v3, "data.frame", info = "fetchVariants with region should return data.frame")
expect_true(nrow(v3) > 0, info = "Region query should return at least one variant")
expect_true(all(v3$CHROM == 22), info = "All variants should be on chr 22")
expect_true(all(v3$POS >= 16050000 & v3$POS <= 16051000),
            info = "All positions should be within region")

# --- getNeighbors() example ---
neighbors <- getNeighbors(ld, "rs587755077", type = "PHASED_R",
                          abs_threshold = sqrt(0.8), genomic_length = 500000)
expect_inherits(neighbors, "character", info = "getNeighbors should return character vector")
expect_true(length(neighbors) >= 1, info = "Should return at least one neighbor")
expect_true("rs587755077" %in% neighbors, info = "Should include the query variant itself")
