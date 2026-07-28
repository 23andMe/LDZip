.check_fetchVariants_inputs <- function(variants) {

  if (!length(variants))
    stop("`variants` must be non-empty")

  if (anyDuplicated(variants))
    stop("`variants` contains duplicate entries")

  if (is.numeric(variants) && any(variants < 1 | variants %% 1 != 0))
    stop("Numeric `variants` indices must be positive integers")

  invisible(TRUE)
}

#' Fetch variant information from an LDZipMatrix object
#'
#' This function retrieves variant metadata including chromosome, position,
#' variant ID, reference allele, and alternate allele from a compressed
#' \code{LDZipMatrix} object. Queries can be made by integer indices (1-based),
#' variant identifiers (IDs), or genomic regions.
#'
#' @param ld An external pointer to an \code{LDZipMatrix} object,
#'   typically created with \code{LDZipMatrix()}.
#' @param variants Integer (1-based), character scalar/vector of variant indices/IDs,
#'   or genomic region string (e.g., "chr1:10000-20000", "22:16050000-16051000").
#'   Genomic regions are inclusive on both start and end positions.
#'
#' @return
#'   A \code{data.frame} with the following columns:
#'   \itemize{
#'     \item \code{idx}: 1-based variant index
#'     \item \code{CHROM}: chromosome
#'     \item \code{POS}: position (base pairs)
#'     \item \code{ID}: variant identifier (e.g., rsID)
#'     \item \code{REF}: reference allele
#'     \item \code{ALT}: alternate allele
#'   }
#'   Results are sorted by variant index.
#'
#' @note
#' \itemize{
#'   \item In order to use variant identifiers or genomic regions (instead of indices),
#'     please index the \code{LDZipMatrix} object using \code{buildIndex(ld)}
#'   \item Variants in the returned data.frame are sorted by variant index
#' }
#'
#' @examples
#' \dontrun{
#' # Open an LDZip matrix
#' prefix <- file.path(system.file("extdata", package = "LDZipMatrix"), "g1k.chr22.ldzip")
#' ld <- LDZipMatrix(prefix)
#'
#' # Build index (if it doesn't exist already)
#' buildIndex(ld)
#'
#' # Query by numeric indices
#' fetchVariants(ld, 1:10)
#'
#' # Query by variant IDs
#' vars <- c("rs587725733", "rs587631919", "rs587661542")
#' fetchVariants(ld, vars)
#'
#' # Query by genomic region
#' fetchVariants(ld, "22:16050000-16051000")
#' }
#'
#' @seealso \code{\link{fetchLD}}, \code{\link{buildIndex}}
#' @export
fetchVariants <- function(ld, variants) {

  if (!inherits(ld, "LDZipMatrix")) {
    stop("`ld` must be an LDZipMatrix object")
  }

  .check_fetchVariants_inputs(variants)

  variant_db_file <- paste(LDZipMatrix_get_prefix_rcpp(ld), "sqlite", sep=".")

  if (!file.exists(variant_db_file)) {
    stop(sprintf("Database file not found: %s\nPlease run buildIndex(ld) first to create the variant index.", variant_db_file))
  }

  # Resolve variant input (indices, rsIDs, or region)
  variant_resolved <- .resolve_variant_input(variants, variant_db_file)
  indices <- variant_resolved$idx

  # Fetch full variant info including CPRA
  result <- get_cpra_from_indices(variant_db_file, indices)

  return(result)
}
