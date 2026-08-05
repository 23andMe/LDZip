#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5) stop("Usage: find_tags_from_vcor.R <vcor> <vars.txt> <variants> <threshold> <output>")

vcor <- read.table(args[1], header = TRUE, stringsAsFactors = FALSE, comment.char = "", check.names = FALSE)
vcor$UNPHASED_R <- as.numeric(vcor$UNPHASED_R)
names(vcor) <- sub("^#", "", names(vcor))

vars <- read.table(args[2], header = FALSE, comment.char = "#", stringsAsFactors = FALSE)
rsid_to_idx <- setNames(0:(nrow(vars)-1), vars[,3])  # column 3 is ID

variant_indices <- as.integer(readLines(args[3]))
threshold <- as.numeric(args[4])

results <- do.call(rbind, lapply(variant_indices, function(var_idx) {
  var_rsid <- vars[var_idx + 1, 3]  # R is 1-based

  idx_a <- which(vcor$ID_A == var_rsid & abs(vcor$UNPHASED_R) >= threshold)
  idx_b <- which(vcor$ID_B == var_rsid & abs(vcor$UNPHASED_R) >= threshold)

  tag_rsids <- c(
    if (length(idx_a) > 0) vcor$ID_B[idx_a] else character(0),
    if (length(idx_b) > 0) vcor$ID_A[idx_b] else character(0)
  )

  tag_indices <- unname(rsid_to_idx[tag_rsids])
  tag_indices <- tag_indices[!is.na(tag_indices) & tag_indices != var_idx]

  if (length(tag_indices) > 0) {
    data.frame(variant = var_idx, tag_variant = tag_indices)
  } else {
    NULL
  }
}))

if (!is.null(results) && nrow(results) > 0) {
  results <- results[order(results$variant, results$tag_variant), ]
  write.table(results, args[5], sep = "\t", row.names = FALSE, quote = FALSE)
} else {
  write("variant\ttag_variant", args[5])
}
