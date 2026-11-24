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
#' @param variants Integer (1-based) or character scalar/vector of variant indices/IDs/CPRAs,
#'   or a genomic region string.
#'   Inputs are resolved in this order:
#'   \itemize{
#'     \item Integer: 1-based variant indices.
#'     \item Character:
#'       \enumerate{
#'         \item A single region string (\code{"chr:start-end"}): all variants in the
#'           region, inclusive of both ends.
#'         \item Otherwise, all are looked up as variant IDs in the ID column.
#'         \item If \emph{none} are found as IDs, all are looked up as CPRA
#'           (\code{"chr:pos:ref:alt"}, e.g. \code{"22:16050115:G:A"}; \code{":"} or
#'           \code{"_"} separators). REF/ALT must match exactly
#'           (alleles are not swapped).
#'       }
#'   }
#' @param allow_missing Logical flag. If \code{TRUE}, variant identifiers that are
#'   not found are dropped with a warning. An error is raised if none are found.
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
fetchVariants <- function(ld, variants, allow_missing = FALSE) {

  if (!inherits(ld, "LDZipMatrix")) {
    stop("`ld` must be an LDZipMatrix object")
  }

  .check_fetchVariants_inputs(variants)

  variant_db_file <- paste(LDZipMatrix_get_prefix_rcpp(ld), "sqlite", sep=".")

  if (!file.exists(variant_db_file)) {
    stop(sprintf("Database file not found: %s\nPlease run buildIndex(ld) first to create the variant index.", variant_db_file))
  }

  # Resolve variant input (indices, rsIDs, or region)
  variant_resolved <- .resolve_variant_input(variants, variant_db_file, allow_missing)
  indices <- variant_resolved$idx

  # Fetch full variant info including CPRA
  result <- get_cpra_from_indices(variant_db_file, indices)

  return(result)
}
