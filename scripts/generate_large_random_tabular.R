generate_random_large_tabular <- function(N,
                                          prefix,
                                          num_nonzero = 100,
                                          min_sep_frac = 0.3) {
  stopifnot(N > 1)
  stopifnot(num_nonzero > 0)
  stopifnot(min_sep_frac > 0 && min_sep_frac < 1)

  min_sep <- ceiling(min_sep_frac * N)

  snp_ids <- paste0("rs", seq_len(N))

  # sample row indices
  i <- sample.int(N, num_nonzero, replace = TRUE)

  # sample offsets guaranteeing far distance
  offset <- sample.int(N - min_sep, num_nonzero, replace = TRUE) + min_sep
  sign <- sample(c(-1L, 1L), num_nonzero, replace = TRUE)

  j <- i + sign * offset

  # wrap into [1, N]
  j <- ((j - 1L) %% N) + 1L

  # enforce upper triangle
  a <- pmin(i, j)
  b <- pmax(i, j)

  x <- runif(num_nonzero, min = -1, max = 1)

    vcor <- data.frame(
    i = a,
    j = b,
    ID_A = snp_ids[a],
    REF_A = "A",
    ALT_A = "C",
    ID_B = snp_ids[b],
    REF_B = "A",
    ALT_B = "C",
    UNPHASED_R = x
    )

    vcor$UNPHASED_R <- formatC(vcor$UNPHASED_R, digits = 6, format = "g")
    vcor <- vcor[order(vcor$i, vcor$j), ]
    vcor$i <- NULL
    vcor$j <- NULL
    colnames(vcor)[1] <- paste0("#", colnames(vcor)[1])

    vars <- data.frame(
        CHROM = 2,
        POS   = seq_len(N) + 10000,
        ID    = snp_ids,
        REF   = "A",
        ALT   = "C"
    )
    
  write.table(vcor,
              file = paste0(prefix, ".vcor"),
              sep = "\t", quote = FALSE, row.names = FALSE)

  write.table(vars,
              file = paste0(prefix, ".vars"),
              sep = "\t", quote = FALSE, row.names = FALSE,
              col.names = c("#CHROM", "POS", "ID", "REF", "ALT"))
}

# CLI entry point
if (sys.nframe() == 0) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) < 2 || length(args) > 4) {
    stop("Usage: Rscript generate_large_random_tabular.R <N> <prefix> [num_nonzero=100] [min_sep_frac=0.3]")
  }

  N <- as.integer(args[1])
  prefix <- args[2]
  num_nonzero <- if (length(args) >= 3) as.integer(args[3]) else 100
  min_sep_frac <- if (length(args) >= 4) as.numeric(args[4]) else 0.3

  generate_random_large_tabular(N, prefix, num_nonzero, min_sep_frac)
}
