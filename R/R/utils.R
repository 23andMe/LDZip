.parse_region <- function(region_str) {
  pattern <- "^(chr)?([0-9XYM]+):([0-9]+)-([0-9]+)$"
  m <- regexec(pattern, region_str)
  matches <- regmatches(region_str, m)[[1]]
  if (length(matches) == 0) stop("Invalid region format. Expected: chr:start-end (e.g., chr1:10000-20000)")
  list(chrom = matches[3], start = as.integer(matches[4]), end = as.integer(matches[5]))
}

.is_region <- function(x) {
  is.character(x) && length(x) == 1 && grepl("^(chr)?[0-9XYM]+:[0-9]+-[0-9]+$", x)
}

.resolve_variant_input <- function(input, variant_db_file) {
  if (.is_region(input)) {
    region <- .parse_region(input)
    r <- get_rsids_by_region(variant_db_file, region)
    idx <- r$idx
    names <- r$rsids
  } else if (is.character(input)) {
    r <- get_indices_by_rsid(variant_db_file, rsids=input)
    idx <- r$idx
    names <- r$rsid
  } else {
    idx <- input
    names <- input
  }
  list(idx = idx, names = names)
}
