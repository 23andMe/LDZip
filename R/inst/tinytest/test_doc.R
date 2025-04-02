library(LDZipMatrix)
library(tinytest)

prefix <- file.path(system.file("extdata", package = "LDZipMatrix"), "g1k.chr22.ldzip")
ld <- LDZipMatrix(prefix)
buildIndex(ld)

# define test sets
idx1  <- 2
idxN  <- 1:5
rsid1 <- "rs587755077"
rsidN <- c("rs587725733","rs587631919","rs587661542","rs587718290","rs587729333")

for (simplify in c(TRUE, FALSE)) {
  for (use_rsid in c(FALSE, TRUE)) {
    for (multi in c(FALSE, TRUE)) {

      one  <- if (use_rsid) rsid1 else idx1
      many <- if (use_rsid) rsidN else idxN
      types <- if (multi) c("PHASED_R", "DPRIME") else "PHASED_R"

      msg <- paste("simplify =", simplify,
                   "use_rsid =", use_rsid,
                   "multi =", multi)

      ## 1x1
      out1 <- fetchLD(ld, one, one, type = types, simplify = simplify)
      if (!multi && simplify) {
        expect_inherits(out1, "numeric", info = msg)
        expect_equal(length(out1), 1, info = msg)
      } else if (multi && simplify) {
        expect_inherits(out1, "data.frame", info = msg)
        expect_equal(nrow(out1), 1, info = msg)
        expect_equal(ncol(out1), length(types), info = msg)
      } else {
        expect_inherits(out1, "list", info = msg)
        expect_equal(length(out1), length(types), info = msg)
        expect_true(all(vapply(out1, is.matrix, logical(1))), info = msg)
      }

      ## 1xN
      out2 <- fetchLD(ld, one, many, type = types, simplify = simplify)
      if (simplify && !multi) {
        expect_true(is.atomic(out2), info = msg)
        expect_equal(length(out2), length(many), info = msg)
      } else if (simplify && multi) {
        expect_inherits(out2, "data.frame", info = msg)
        expect_equal(nrow(out2), length(many), info = msg)
        expect_equal(ncol(out2), length(types), info = msg) # row+col ids + stats
      } else {
        expect_inherits(out2, "list", info = msg)
        expect_equal(length(out2), length(types), info = msg)
        expect_true(all(vapply(out2, is.matrix, logical(1))), info = msg)
      }

      ## NxN
      out3 <- fetchLD(ld, many, many, type = types, simplify = simplify)
      if (simplify && !multi) {
        expect_true(is.matrix(out3), info = msg)
        expect_equal(dim(out3), c(length(many), length(many)), info = msg)
      } else if (simplify && multi) {
        expect_inherits(out3, "list", info = msg)
        expect_equal(length(out3), length(types), info = msg)
        expect_true(all(vapply(out3, is.matrix, logical(1))), info = msg)
      } else {
        expect_inherits(out3, "list", info = msg)
        expect_equal(length(out3), length(types), info = msg)
        expect_true(all(vapply(out3, is.matrix, logical(1))), info = msg)
      }
    }
  }
}


# --- getNeighbors sanity checks ---------------------------------------------

# RSID input → RSID output
res1 <- getNeighbors(ld, "rs587755077",
                     type = "PHASED_R",
                     abs_threshold = sqrt(0.8),
                     genomic_length = 500000)
expect_inherits(res1, "character")
expect_true(length(res1) >= 1)

# Index input → Index output
res_idx <- getNeighbors(ld, 2,
                        type = "PHASED_R",
                        abs_threshold = sqrt(0.8),
                        genomic_length = 500000)
expect_true(is.integer(res_idx))
expect_true(length(res_idx) >= 1)

# RSID input for D' perfect LD
res2 <- getNeighbors(ld, "rs587755077",
                     type = "DPRIME",
                     abs_threshold = 1,
                     genomic_length = 100000)
expect_inherits(res2, "character")

# Verify D' values are 1
out <- fetchLD(ld, "rs587755077", res2, type = "DPRIME", simplify = TRUE)
if (length(out) > 0) {
  if (is.data.frame(out)) {
    expect_true(all(abs(out$DPRIME - 1) < 1e-8))
  } else if (is.numeric(out)) {
    expect_true(all(abs(out - 1) < 1e-8))
  }
}

# --- Pairwise queries -------------------------------------------------------

# Pairwise with integer indices (single statistic)
rows <- c(1, 5, 10)
cols <- c(2, 8, 15)
pw1 <- fetchLD(ld, rows, cols, types = "PHASED_R", pairwise = TRUE)
expect_inherits(pw1, "data.frame", info = "Pairwise should return data.frame")
expect_equal(nrow(pw1), 3, info = "Pairwise should have 3 rows")
expect_equal(ncol(pw1), 3, info = "Single stat pairwise should have 3 columns (stat + var1 + var2)")
expect_true("PHASED_R" %in% colnames(pw1), info = "Column should be named PHASED_R")
expect_true(all(c("var1", "var2") %in% colnames(pw1)), info = "Should have var1 and var2 columns")

# Pairwise with variant IDs
row_vars <- c("rs587725733", "rs587631919", "rs587661542")
col_vars <- c("rs587755077", "rs587718290", "rs587729333")
pw2 <- fetchLD(ld, row_vars, col_vars, types = "PHASED_R", pairwise = TRUE)
expect_inherits(pw2, "data.frame", info = "Pairwise with rsids should return data.frame")
expect_equal(nrow(pw2), 3, info = "Pairwise with rsids should have 3 rows")
expect_equal(ncol(pw2), 3, info = "Single stat pairwise should have 3 columns (stat + var1 + var2)")

# Pairwise with multiple statistics
pw3 <- fetchLD(ld, rows, cols, types = c("PHASED_R", "DPRIME"), pairwise = TRUE)
expect_inherits(pw3, "data.frame", info = "Multi-stat pairwise should return data.frame")
expect_equal(nrow(pw3), 3, info = "Multi-stat pairwise should have 3 rows")
expect_equal(ncol(pw3), 4, info = "Multi-stat pairwise should have 4 columns (2 stats + var1 + var2)")
expect_true(all(c("PHASED_R", "DPRIME") %in% colnames(pw3)),
            info = "Columns should be named PHASED_R and DPRIME")

# Verify pairwise values match non-pairwise queries and are in input order
for (i in seq_along(rows)) {
  r <- rows[i]
  c <- cols[i]
  val_regular <- fetchLD(ld, r, c, types = "PHASED_R", simplify = TRUE)
  # Pairwise returns results in same order as input
  expect_equal(pw1$PHASED_R[i], val_regular, tolerance = 1e-10,
               info = sprintf("Pairwise value for (%d,%d) should match regular query", r, c))
  expect_equal(pw1$var1[i], r, info = sprintf("var1[%d] should be %d", i, r))
  expect_equal(pw1$var2[i], c, info = sprintf("var2[%d] should be %d", i, c))
}

# Test duplicate pairs are allowed and return correct values
rows_dup <- c(1, 1, 5)
cols_dup <- c(2, 2, 8)
pw_dup <- fetchLD(ld, rows_dup, cols_dup, types = "PHASED_R", pairwise = TRUE)
expect_equal(nrow(pw_dup), 3, info = "Should return all pairs including duplicates")
expect_equal(pw_dup$PHASED_R[1], pw_dup$PHASED_R[2], info = "Duplicate pairs should have same value")
expect_equal(pw_dup$var1[1], pw_dup$var1[2], info = "Duplicate pairs should have same var1")
expect_equal(pw_dup$var2[1], pw_dup$var2[2], info = "Duplicate pairs should have same var2")

# Test N×1 query (in addition to 1×N)
out_Nx1 <- fetchLD(ld, idxN, idx1, types = "PHASED_R", simplify = TRUE)
expect_true(is.atomic(out_Nx1), info = "N×1 single stat should return vector")
expect_equal(length(out_Nx1), length(idxN), info = "N×1 should have length N")

# Test that non-pairwise results are sorted by indices (not input order)
unsorted_idx <- c(5, 1, 3, 2, 4)
out_unsorted <- fetchLD(ld, unsorted_idx, unsorted_idx, types = "PHASED_R", simplify = TRUE)
expect_true(is.matrix(out_unsorted), info = "Should return matrix")
expect_equal(rownames(out_unsorted), as.character(sort(unsorted_idx)), info = "Row names should be sorted")
expect_equal(colnames(out_unsorted), as.character(sort(unsorted_idx)), info = "Column names should be sorted")

