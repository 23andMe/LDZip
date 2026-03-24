library(LDZipMatrix)
library(tinytest)

source("../../../scripts/source.R")

N <- 1000

test_cases <- list(
  list(name="single_snp",   rsid="rs20",   window=1000),
  list(name="edge_snp",     rsid="rs1",    window=500),
  list(name="end_snp",     rsid="rs1000",    window=500),
  list(name="middle_snp",   rsid="rs500",  window=2000),
  list(name="dense_region", rsid="rs100",  window=1000)
)

# load original matrix
df <- read.table("../../tests/data/test_tabular.vcor", header = T, stringsAsFactors = FALSE)


# list of compressed outputs
bits_list <- c(8, 16, 32, 99)
min_list  <- c(0.00, 0.001, 0.01, 0.1)

for (bits in bits_list) {
  for (min in min_list) {
    tol <- get_threshold(bits, min)

    prefix <- sprintf("../../tests/data/compressed_bits_tabular_%d_min_%s", bits, min)
    cat("Testing:", prefix, "\n")
    
    # load compressed
    ld <- LDZipMatrix(prefix)
    suppressMessages(buildIndex(ld))

    for (tc in test_cases) {
        rsid=tc$rsid

        min_thresh = 0.0
        df_subset = subset(df, abs(UNPHASED_R) > as.numeric(min) & (df$UNPHASED_R^2) >=min_thresh)
        val_old <- unique(c(df_subset[df_subset$ID_A == rsid, "ID_B"], 
                          df_subset[df_subset$ID_B == rsid, "ID_A"], 
                          rsid))   # add diagonal entry
        val_old <- sort(val_old)

        ## load new data
        val_new <- sort(getNeighbors(ld, rsid, type="UNPHASED_R", abs_threshold=sqrt(min_thresh), genomic_length=100000))

        expect_true(
          setequal(val_new, val_old),
          info = sprintf("Mismatch for prefix %s", prefix)
        )

    }

  }
}

