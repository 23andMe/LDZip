library(LDZipMatrix)
library(tinytest)

source("../../../scripts/source.R")

N <- 1000

test_cases <- list(
    # one to one
    list(name="first_first",        type="one_to_one",  category="diagonal",      row=1,          col=1),
    list(name="last_last",          type="one_to_one",  category="diagonal",      row=N,          col=N),
    list(name="off_diag_small",     type="one_to_one",  category="off-diagonal",  row=10,         col=20),
    list(name="off_diag_large",     type="one_to_one",  category="off-diagonal",  row=900,        col=600),
    list(name="first_last",         type="one_to_one",  category="edge",          row=1,          col=N),
    list(name="last_first",         type="one_to_one",  category="edge",          row=N,          col=1),
    list(name="last_near_diag",     type="one_to_one",  category="off-diagonal",  row=N,          col=N-1),

    # one to many
    list(name="range",              type="one_to_many", category="row-fixed",     row=1,          col=1:100),
    list(name="range_rev",          type="one_to_many", category="col-fixed",     row=1:100,      col=50),
    list(name="random_row_subset",  type="one_to_many", category="row-random",    row=sort(sample(1:N, N/5)), col=50),
    list(name="random_col_subset",  type="many_to_one", category="col-random",    row=50,         col=sort(sample(1:N, N/5))),

    # many to many
    list(name="upper_triangle_slice", type="many_to_many", category="upper-tri",   row=100:120,  col=800:820),
    list(name="wide_row_band",        type="many_to_many", category="row-band",    row=200:210,  col=c(10,20,40,50)),
    list(name="full_quadrant",        type="many_to_many", category="block",       row=1:100,    col=901:1000),
    list(name="rand_block_square",    type="many_to_many", category="block",       row=sort(sample(1:N, 20)),  col=sort(sample(1:N, 20))),
    list(name="rand_row_band",        type="many_to_many", category="row-band",    row=sort(sample(1:N, 15)),  col=sort(sample(1:N, 5))),
    list(name="rand_col_band",        type="many_to_many", category="col-band",    row=sort(sample(1:N, 5)),   col=sort(sample(1:N, 15))),
    list(name="rand_upper_triangle",  type="many_to_many", category="upper-tri",   row=sort(sample(1:(N/2), 10)), col=sort(sample((N/2):N, 10))),
    list(name="rand_lower_triangle",  type="many_to_many", category="lower-tri",   row=sort(sample((N/2):N, 10)), col=sort(sample(1:(N/2), 10)))
)
# load original matrix
orig <- read_tabular_ld_matrix("../../tests/data/test_tabular")

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
        r <- paste("rs",tc$row,sep="") 
        c <- paste("rs",tc$col,sep="")

        val_new <- fetchLD(ld, r, c, type="UNPHASED_R", simplify=F)$UNPHASED_R
        val_old <- orig[r, c]

        expect_equal(
            as.vector(val_new),
            as.vector(val_old),
            tolerance = tol + 1e-7,
            scale = 1,
            info = sprintf("Mismatch for prefix %s", prefix)
        )
    }

  }
}
