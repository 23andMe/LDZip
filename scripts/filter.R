#!/usr/bin/env Rscript

library(argparse)
library(jsonlite)
source("../../scripts/source.R")

# Argument parser
args=list(orig = "output/plink_tabular.vcor", keep="output/random_keep", out="output/bla")

parser <- ArgumentParser(description = "Filter LD data by index list")
parser$add_argument("-orig", type = "character", required = TRUE,
                    help = "Original LD file (.bin or .vcor)")
parser$add_argument("-keep", type = "character", required = TRUE,
                    help = "File with 0-based indices to keep")
parser$add_argument("-out", type = "character", required = TRUE,
                    help = "Output LD file")

args <- parser$parse_args()

ld_type <- detect_ld_type(args$orig)

# Read keep indices (0-based → 1-based)
keep_idx <- scan(args$keep, what = integer(), quiet = TRUE) + 1

# -------- VCOR path --------
# -------- VCOR path --------
if (ld_type == "vcor") {

    vars_file <- paste0(args$orig, ".vars.txt")
    if (!file.exists(vars_file))
        stop("Missing variant file: ", vars_file)

vars <- read.table(
    vars_file,
    header = TRUE,
    comment.char = "",
    stringsAsFactors = FALSE,
    check.names = FALSE
)
    # Validate keep indices
    keep_idx <- keep_idx[keep_idx >= 1 & keep_idx <= nrow(vars)]
    if (length(keep_idx) == 0)
        stop("No valid keep indices after bounds check")
print(nrow(vars))
# print(keep_idx)
print(colnames(vars))
    keep_ids <- vars$ID[keep_idx]
    keep_ids <- keep_ids[!is.na(keep_ids)]

    if (length(keep_ids) == 0)
        stop("All keep IDs are NA")

    # Read header
    vcor_header <- readLines(args$orig, n = 1)

    # Read VCOR data (skip header only)
    vcor <- read.table(
        args$orig,
        header = FALSE,
        skip = 1,
        sep = "",
        stringsAsFactors = FALSE
    )

    colnames(vcor) <- strsplit(sub("^#", "", vcor_header), "\\s+")[[1]]

    # Subset rows where BOTH variants are kept
    vcor_sub <- vcor[
        vcor$ID_A %in% keep_ids & vcor$ID_B %in% keep_ids,
        ,
        drop = FALSE
    ]

    if (nrow(vcor_sub) == 0)
        stop("VCOR filtering produced 0 rows — check keep file vs vars")

    # Write header
    cat(vcor_header, file = args$out, sep = "\n")

    # Write rows
    write.table(
        vcor_sub,
        args$out,
        quote = FALSE,
        row.names = FALSE,
        col.names = FALSE,
        append = TRUE,
        sep="\t"
    )

    # Write subset vars
    write.table(
        vars[keep_idx, , drop = FALSE],
        paste0(args$out, ".vars.txt"),
        quote = FALSE,
        row.names = FALSE,
        col.names = TRUE
    )

    cat("Filtered VCOR written to", args$out, "\n")
    quit(status = 0)
}


# -------- Binary matrix path --------

# Read matrix
old_val <- read_ld_data(args$orig)
N <- sqrt(length(old_val))
stopifnot(N == as.integer(N))
old_mat <- matrix(old_val, nrow = N, ncol = N, byrow = FALSE)

# Subset
new_mat <- old_mat[keep_idx, keep_idx, drop = FALSE]

# Write binary output
writeBin(as.numeric(new_mat), args$out, size = 4)

# ---- handle .vars file if present ----
vars_file <- paste0(args$orig, ".vars")
if (file.exists(vars_file)) {
    vars <- read.table(vars_file, header = FALSE, stringsAsFactors = FALSE)
    vars_sub <- vars[keep_idx, , drop = FALSE]
    write.table(
        vars_sub,
        paste0(args$out, ".vars"),
        quote = FALSE,
        row.names = FALSE,
        col.names = FALSE
    )
}

cat("Filtered LD matrix written to", args$out, "\n")
