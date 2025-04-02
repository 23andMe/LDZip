#' Find neighboring variants in an LDZipMatrix
#'
#' This function finds all variants exceeding a specified LD threshold
#' for a given statistics type, within a given genomic window around 
#' a query variant, using an \code{LDZipMatrix} object.
#'
#' @param ld An external pointer to an \code{LDZipMatrix} object,
#'   typically created with \code{LDZipMatrix()}.
#' @param variant Integer (1-based)/Character scalar referring to the reference variant index/identifier around
#'   which neighbors are sought.
#' @param type Character scalar referring to the type of linkage statistic sought.
#' @param genomic_length Integer scalar referring to the window size in base pairs
#'   (e.g., 1e6 for +/-1 Mb).
#'
#' @return Integer/Character vector of variant indices/indentifiers. 
#' \itemize{
#' \item If \code{variant} is a variant index (\code{Integer}), the function returns an integer vector (1-based) of variant indices
#' \item If \code{variant} is a variant identifier (\code{Character}), the function returns a character vector of variant identifiers
#' }
#' @note
#' \itemize{
#'   \item In order to use variant identifiers (instead of indices), please index the \code{LDZipMatrix} object using \code{buildIndex(ld)}
#' }
#' 
#' @examples
#' \dontrun{
#' 
#' # Open an LDZip matrix
#' prefix <- file.path(system.file("extdata", package = "LDZipMatrix"), "g1k.chr22.ldzip")
#' ld <- LDZipMatrix(prefix)
#' 
#' # Build index (if it doesn't exist already)
#' buildIndex(ld)
#' 
#' # Find neighbors with r^2 > 0.8 within +/-500kb of rs587755077
#' getNeighbors(ld, "rs587755077", type="PHASED_R", abs_threshold = sqrt(0.8), genomic_length=500000)
#' 
#' 
#' # Find neighbors within +/-100Kb of rs587755077 that are in perfect LD (D'=1)
#' getNeighbors(ld, "rs587755077", type="DPRIME", abs_threshold = 1, genomic_length=100000)
#' 
#' # Verify their D' value
#' res = getNeighbors(ld, "rs587755077", type="DPRIME", abs_threshold = 1.0, genomic_length=100000)
#' fetchLD(ld, "rs587755077", res, "DPRIME") 
#' 
#' }
#'
#' @export
getNeighbors <- function(ld, variant, type="UNPHASED_R2", abs_threshold=0.8, genomic_length=1e5) {
  variant_db_file <- paste(LDZipMatrix_get_prefix_rcpp(ld), "sqlite", sep=".")
  if (!file.exists(variant_db_file)) {
    stop(sprintf("Database file not found: %s", variant_db_file))
  }
  if (is.character(variant)) {
    indices <- get_indices_by_rsid(variant_db_file, rsids = variant)
    variant_index  <- indices$idx
    cpra <- get_cpra_from_rsids(variant_db_file, variant)
  }else {
    variant_index = variant
    cpra <- get_cpra_from_indices(variant_db_file, variant_index)
  }

  if (nrow(cpra) == 0) stop("variant not found in database")

  chrom <- cpra$CHROM[1]; pos <- cpra$POS[1]
  region <- list(chrom = chrom, start = pos - genomic_length, end = pos + genomic_length)
  nnz_idx <- LDZipMatrix_get_neighbors_rcpp(ld, variant_index-1L, abs_threshold, type) + 1L

  if (is.character(variant)) {
    variants_in_region <- get_rsids_by_region(variant_db_file, region)
    idx_intersect <- which(variants_in_region$idx %in%  nnz_idx)
    return (variants_in_region$rsids[idx_intersect])
  } else {
    variants_in_region <- get_indices_by_region(variant_db_file, region)
    idx_intersect <- intersect(variants_in_region, nnz_idx)
    return (idx_intersect)
  }

}