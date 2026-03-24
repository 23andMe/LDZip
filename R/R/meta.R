#' @name LDZipMatrix-metadata
#' @title Functions for LDZipMatrix metadata
#' @aliases nrow.LDZipMatrix ncol.LDZipMatrix dim.LDZipMatrix nnz stats bits layout
#' 
#' @description 
#' Retrieve dimensions and metadata from an LDZipMatrix object
#'
#' These functions provide convenient accessors for structural information
#' and metadata associated with a compressed \code{LDZipMatrix} object.
#' They mirror base R generics (e.g. \code{nrow()}, \code{ncol()}, \code{dim()})
#' where appropriate, and add specialized accessors for sparse matrix properties
#' and LD-specific metadata.
#'
#' @param x An external pointer to an \code{LDZipMatrix} object,
#'   typically created with \code{LDZipMatrix()}.
#' 
#' @details
#' \itemize{
#'   \item \code{nrow(x)}: Return the number of rows (variants) in the matrix.
#'   \item \code{ncol(x)}: Return the number of columns (variants) in the matrix.
#'   \item \code{dim(x)}: Return both number of rows and columns as a vector.
#'   \item \code{nnz(x)}: Return the number of non-zero entries stored in the compressed representation.
#'   \item \code{stats(x)}: Return the set of LD statistics available in the object.
#'     Possible values include \code{"PHASED_R"}, \code{"UNPHASED_R"},
#'     \code{"PHASED_R2"}, \code{"UNPHASED_R2"}, \code{"D"}, \code{"DPRIME"}.
#'   \item \code{bits(x)}: Return the number of bits used for compression (e.g. 8, 16, 32).
#'   \item \code{layout(x)}: Return the matrix layout used in compression (e.g. \code{"upper"}, \code{"full"}).
#' }
#'
#' These accessors are lightweight and do not decompress the full LD matrix.
#' They are intended for quick inspection and for programmatic use within
#' higher-level analysis workflows.
#'
#' @return
#' An integer, logical, or character vector depending on the function:
#' \itemize{
#'   \item \code{nrow()}, \code{ncol()} return integers.
#'   \item \code{dim()} returns an integer vector of length two.
#'   \item \code{nnz()} returns an integer.
#'   \item \code{stats()} returns a character vector of available LD statistics.
#'   \item \code{bits()} returns an integer.
#'   \item \code{layout()} returns a character string.
#' }
#'
#' @examples
#' \dontrun{
#' prefix <- file.path(system.file("extdata", package = "LDZipMatrix"), "g1k.chr22.ldzip")
#'
#' # Create a new LDZipMatrix handle
#' ld <- LDZipMatrix(prefix)
#' 
#' # Query metadata
#' nrow(ld)
#' ncol(ld)
#' dim(ld)
#' nnz(ld)
#' stats(ld)
#' bits(ld)
#' layout(ld)
#' print(ld)
#' summary(ld)
#' }
#'
#' @seealso \code{\link{LDZipMatrix}}

#' @rdname LDZipMatrix-metadata
#' @export
#' @method nrow LDZipMatrix
nrow.LDZipMatrix <- function(x) {
  LDZipMatrix_nrows_rcpp(x)
}

#' @rdname LDZipMatrix-metadata
#' @export
#' @method ncol LDZipMatrix
ncol.LDZipMatrix <- function(x) {
  LDZipMatrix_ncols_rcpp(x)
}


#' @rdname LDZipMatrix-metadata
#' @export
dim.LDZipMatrix <- function(x) {
  c(LDZipMatrix_nrows_rcpp(x), LDZipMatrix_ncols_rcpp(x))
}


#' @rdname LDZipMatrix-metadata
#' @export
print.LDZipMatrix <- function(x, ...) {
  d <- dim(x)
  cat("<LDZipMatrix>\n")
  cat("  Dimensions:", d[1], "x", d[2], "\n")
  cat("  Stats     :", paste(stats(x), collapse = ", "), "\n")
  invisible(x)
}


#' @rdname LDZipMatrix-metadata
#' @export
summary.LDZipMatrix <- function(object, ...) {
  d <- dim(object)
  out <- list(
    dimensions = d,
    nnz        = nnz(object),
    bits       = bits(object),
    layout     = LDZipMatrix_format_rcpp(object),
    stats      = stats(object)
  )
  class(out) <- "summary.LDZipMatrix"
  out
}


#' @rdname LDZipMatrix-metadata
#' @export
print.summary.LDZipMatrix <- function(x, ...) {
  cat("<Summary of LDZipMatrix>\n")
  cat("  Dimensions:", x$dimensions[1], "x", x$dimensions[2], "\n")
  cat("  Nonzeros  :", x$nnz, "\n")
  cat("  Bits      :", x$bits, "\n")
  cat("  Layout    :", x$layout, "\n")
  cat("  Stats     :", paste(x$stats, collapse = ", "), "\n")
  invisible(x)
}

#' @rdname LDZipMatrix-metadata
#' @export
nnz <- function(x) {
  LDZipMatrix_nnz_rcpp(x)
}

#' @rdname LDZipMatrix-metadata
#' @export
layout <- function(x, ...) {
  LDZipMatrix_format_rcpp(x)
}


#' @rdname LDZipMatrix-metadata
#' @export
stats <- function(x) {
  all_stats <- LDZipMatrix_allStats_rcpp()
  all_stats[
    vapply(all_stats, function(si) LDZipMatrix_hasStatString_rcpp(x, si), logical(1))
  ]
}


#' @rdname LDZipMatrix-metadata
#' @export
bits <- function(x) {
  LDZipMatrix_bits_rcpp(x)
}
