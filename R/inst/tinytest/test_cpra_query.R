library(LDZipMatrix)
library(tinytest)

# CPRA (chr:pos:ref:alt) fallback when no identifier is found in the ID column

ld <- LDZipMatrix(file.path(system.file("extdata", package = "LDZipMatrix"), "g1k.chr22.ldzip"))
suppressMessages(buildIndex(ld))

rsids <- c("rs587755077", "rs587631919", "rs587661542", "rs587718290")
cpras <- c("22:16050115:G:A", "22:16058342:A:T", "22:16059796:A:C", "22:16061912:C:A")
q <- function(expr) suppressWarnings(expr)

# CPRA gives the same result as rsID
expect_warning(fetchVariants(ld, cpras), "matched by CPRA")
expect_equal(q(fetchVariants(ld, cpras)), fetchVariants(ld, rsids))
expect_equal(unname(q(fetchLD(ld, cpras[1], cpras[2:4], "PHASED_R"))), unname(fetchLD(ld, rsids[1], rsids[2:4], "PHASED_R")))
expect_equal(q(fetchLD(ld, cpras[1:3], cpras[2:4], "PHASED_R", pairwise = TRUE))$PHASED_R,
             fetchLD(ld, rsids[1:3], rsids[2:4], "PHASED_R", pairwise = TRUE)$PHASED_R)

# accepted formats: chr prefix, '_' separator, lower-case alleles, sex chromosome
expect_equal(q(fetchVariants(ld, c("chr22:16050115:G:A", "22_16058342_a_t")))$ID, rsids[1:2])
ld_y <- LDZipMatrix("../../tests/data/compressed_chrY")
suppressMessages(buildIndex(ld_y))
expect_equal(q(fetchVariants(ld_y, "Y:10000:A:C"))$ID, "rs1")

# errors
expect_error(fetchVariants(ld, "22:16050115:A:G"), "No matching variants")   # alleles not swapped
expect_error(fetchVariants(ld, c(cpras[1], "22:1:A:G")), "Missing CPRAs")
expect_error(fetchVariants(ld, c(rsids[1], cpras[2])), "Missing rsIDs")       # no CPRA fallback if some IDs match
expect_error(fetchVariants(ld, "rs99999999"), "Could not parse")

# allow_missing: drops some, never all
expect_equal(q(fetchVariants(ld, c(rsids[1:2], "rs99999999"), allow_missing = TRUE))$ID, rsids[1:2])
expect_equal(q(fetchVariants(ld, c(cpras[1:2], "22:1:A:G"), allow_missing = TRUE))$ID, rsids[1:2])
expect_equal(q(fetchLD(ld, rsids[1], c(rsids[2:4], "rs99999999"), "PHASED_R", allow_missing = TRUE)),
             fetchLD(ld, rsids[1], rsids[2:4], "PHASED_R"))
expect_error(fetchVariants(ld, "22:1:A:G", allow_missing = TRUE), "No matching variants")
expect_error(fetchLD(ld, rsids[1:2], rsids[3:4], pairwise = TRUE, allow_missing = TRUE), "not supported")
