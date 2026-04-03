generate_random_tabular <- function(N, prefix, type, num_nonzero = 10000, chrom = "2") {
  # start with zeros + identity
  mat <- matrix(0, N, N)
  diag(mat) <- 1

  # pick a few random off-diagonal pairs
  idx <- which(upper.tri(mat), arr.ind = TRUE)
  chosen <- idx[sample(nrow(idx), min(num_nonzero, nrow(idx))), , drop = FALSE]
  for (k in seq_len(nrow(chosen))) {
    i <- chosen[k, 1]
    j <- chosen[k, 2]
    val <- runif(1, min = -1, max = 1)
    mat[i, j] <- val
    mat[j, i] <- val
  }

  # SNP info
  snp_ids <- paste0("rs", 1:N)
  pos <- seq(10000, 10000 + N - 1)
  ref <- rep("A", N)
  alt <- rep("C", N)

  vars <- data.frame(CHROM = chrom, POS = pos, ID = snp_ids, REF = ref, ALT = alt)

  # pairwise values
  rows <- list()
  for (i in 1:(N - 1)) {
    for (j in (i + 1):N) {
      val <- mat[i, j]
      if (val != 0) {
        rows[[length(rows) + 1]] <- data.frame(
          ID_A = snp_ids[i],
          REF_A = "A",
          ALT_A = "C",
          ID_B = snp_ids[j],
          REF_B = "A",
          ALT_B = "C",
          setNames(list(val), type)
        )
      }
    }
  }

  vcor <- do.call(rbind, rows)

  # write files
  vcor_names <- names(vcor)
  vcor_names[1] <- paste0("#", vcor_names[1])
  write.table(vcor,
              file = paste0(prefix, ".vcor"),
              sep = "\t", quote = FALSE, row.names = FALSE,
              col.names = vcor_names)
  write.table(vars,
              file = paste0(prefix, ".vars"),
              sep = "\t", quote = FALSE, row.names = FALSE,
              col.names = c("#CHROM", "POS", "ID", "REF", "ALT"))

  invisible(list(mat = mat, vars = vars, vcor = vcor))
}

# only run this block when called via `Rscript`
if (sys.nframe() == 0) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) < 2 || length(args) > 4) {
    stop("Usage: Rscript generate_random_tabular.R <N> <output_prefix> [type] [chrom]")
  }
  type <- if (length(args) >= 3) args[3] else "UNPHASED_R"
  chrom <- if (length(args) >= 4) args[4] else "2"
  generate_random_tabular(as.integer(args[1]), args[2], type, chrom = chrom)
}
