library(LDZipMatrix)
library(tinytest)

source("../../../scripts/source.R")

N <- 1000

# Test that chr Y data works correctly with getNeighbors
test_cases <- list(
  list(name="single_snp_chrY",   rsid="rs20",   window=1000),
  list(name="edge_snp_chrY",     rsid="rs1",    window=500),
  list(name="middle_snp_chrY",   rsid="rs500",  window=2000)
)

# Load original chr Y matrix
df <- read.table("../../tests/data/test_chrY.vcor", header = T, stringsAsFactors = FALSE, comment.char = "", check.names = FALSE)

# Compressed output for chr Y
prefix <- "../../tests/data/compressed_chrY"
cat("Testing chr Y:", prefix, "\n")

# Load compressed
ld <- LDZipMatrix(prefix)
suppressMessages(buildIndex(ld))

# Verify chromosome is stored as Y
vars <- read.table(paste0(prefix, ".vars.txt"), header = TRUE, stringsAsFactors = FALSE)
expect_true(
  all(vars$CHROM == "Y"),
  info = "Chromosome should be Y"
)

for (tc in test_cases) {
    rsid <- tc$rsid

    min_thresh <- 0.0
    df_subset <- subset(df, abs(UNPHASED_R) > 0.0 & (df$UNPHASED_R^2) >= min_thresh)
    val_old <- unique(c(df_subset[df_subset$`#ID_A` == rsid, "ID_B"],
                      df_subset[df_subset$ID_B == rsid, "#ID_A"],
                      rsid))   # add diagonal entry
    val_old <- sort(val_old)

    ## load new data using getNeighbors
    val_new <- sort(getNeighbors(ld, rsid, type="UNPHASED_R", abs_threshold=sqrt(min_thresh), genomic_length=100000))

    expect_true(
      setequal(val_new, val_old),
      info = sprintf("Mismatch for rsid %s on chr Y", rsid)
    )
}

cat("✓ Chr Y test passed\n")
