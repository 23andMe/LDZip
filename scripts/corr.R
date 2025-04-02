#!/usr/bin/env Rscript

library(argparse)
library(jsonlite)
source("../../scripts/source.R")

args=list(orig = c("output/chr20.unphased.vcor1.bin", "output/chr21.unphased.vcor1.bin", "output/chr22.unphased.vcor1.bin"), new="output/concat_full.bin")
args=list(orig = c("output/chr20.vcor", "output/chr21.vcor", "output/chr22.vcor"), new="output/concat_full.vcor")
args=list(orig = c("output/test1.bin", "output/test2.bin", "output/test3.bin", "output/test4.bin"), new="output/concat_full.bin")

# Argument parser
parser <- ArgumentParser(description = "Compare roundtrip LD matrices")
parser$add_argument("-orig", type = "character", nargs = "+", required = TRUE,
                    help = "Original LD file(s), space separated")
parser$add_argument("-new", type = "character", nargs = "+", required = TRUE,
                    help = "New LD file(s), space separated")
# parser$add_argument("-col", type = "character", nargs = "+", required = TRUE,
#                     help = "New LD file(s), space separated")
args <- parser$parse_args()


if(detect_ld_type(args$orig) == "vcor" && detect_ld_type(args$new) == "vcor") compare_ld_headers(args$orig, args$new)



old_val <- read_ld_data(args$orig)
new_val <- read_ld_data(args$new)
stopifnot(length(old_val) == length(new_val))

# Metrics
num_non_zero_old <- sum(old_val != 0.0, na.rm=TRUE)
num_non_zero_new <- sum(new_val != 0.0, na.rm=TRUE)
cor_full         <- cor(old_val, new_val, use = "complete.obs")
cor_nonzero      <- cor(old_val[old_val != 0], new_val[old_val != 0], use = "complete.obs")
max_diff         <- max(abs(old_val - new_val), na.rm=TRUE)

stats <- list(
  num_non_zero_old = num_non_zero_old,
  num_non_zero_new = num_non_zero_new,
  cor_full = cor_full,
  cor_nonzero = cor_nonzero,
  max_diff = max_diff
)

write_json(stats, paste0( tools::file_path_sans_ext(args$new), "_stats.json"), pretty = TRUE, auto_unbox = TRUE, digits = 6)
cat("Stats written to ", args$new, "_stats.json\n")
