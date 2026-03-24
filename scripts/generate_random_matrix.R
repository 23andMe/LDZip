generate_random_matrix <- function(N, prefix) {
  mat <- matrix(0, N, N)
  mat[upper.tri(mat)] <- runif(length(mat[upper.tri(mat)]), min = -1, max = 1)
  mat[lower.tri(mat)] <- t(mat)[lower.tri(mat)]
  diag(mat) <- 1

  writeBin(as.vector(mat), paste0(prefix, ".bin"), size = 4)
  write.table(1:N,
              paste0(prefix, ".bin.vars"),
              quote = FALSE, row.names = FALSE, col.names = FALSE)
  invisible(mat)
}

# only run this block when called via `Rscript`
if (sys.nframe() == 0) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) != 2) {
    stop("Usage: Rscript gen_matrix.R <N> <output_prefix>")
  }
  generate_random_matrix(as.integer(args[1]), args[2])
}
