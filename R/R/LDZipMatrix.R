#' @name LDZipMatrix-class
#' @title LDZipMatrix Class and Constructor
#'
#' @description
#' An \code{LDZipMatrix} stores pairwise linkage disequilibrium (LD) statistics
#' between genetic variants using a compressed representation on disk.  A single
#' \code{LDZipMatrix} file can also hold multiple LD statistics (e.g., \eqn{r}, \eqn{r^2},
#' \eqn{D}, \eqn{D'}). Instead of recalculating LD on demand, results can be stored once 
#' and retrieved whenever necessary. Because LD matrices grow quadratically with 
#' the number of variants, they are prohibitively large to hold in dense form. 
#' \code{LDZipMatrix} stores them in compressed format that minimizes storage 
#' while still supporting efficient random access to rows, columns, or individual elements as needed.

#' @details
#' The \code{LDZipMatrix} class is implemented as an \emph{external pointer} to a
#' C++ object. The underlying C++ object manages access to a set of compressed
#' files on disk that together encode the matrix in \emph{Compressed Sparse Column (CSC)}
#' format. This format stores, for each column, a pointer to its nonzero entries
#' (\code{.p.bin}), the corresponding row indices (\code{.i.bin}), and one or more
#' compressed value streams (\code{.x.<stats>.bin}). In fact, it stores only the
#' delta encoded values of row indices in (\code{.i.bin}) so that they fit within 
#' a 16 bit integer, and stores the few overflowing delta values in (\code{.io.bin}) 
#'
#' Each \code{.x.<stats>.bin} file corresponds to a particular LD statistic
#' (e.g., \code{PHASED_R}, \code{UNPHASED_R2}, \code{D}, \code{DPRIME}). Multiple
#' such files may coexist for the same matrix prefix, allowing the same structural
#' representation (\code{i}, \code{p}) to be reused across different LD formats.
#'
#' The full LD matrix is never loaded into memory; instead, only the requested rows,
#' columns, or elements are decompressed on demand. This design enables efficient
#' random access to very large LD matrices without requiring gigabytes of RAM.

#'
#' @note
#' A valid LDZipMatrix requires several files with the same \code{prefix}:
#' \itemize{
#'   \item \code{<prefix>.i.bin} – row indices for nonzero entries delta encoded
#'   \item \code{<prefix>.io.bin} – overflowing row indices for nonzero entries
#'   \item \code{<prefix>.p.bin} – column pointers marking the start of each column in \code{i} and \code{x}
#'   \item \code{<prefix>.x.<stats>.bin} – compressed numeric values for a given LD statistic
#'   \item \code{<prefix>.meta.json} – json file describing associated metadata
#'   \item \code{<prefix>.vars.txt} – variant identifiers
#'   \item (optionally) \code{<prefix>.sqlite} – SQLite database with indexed variant information
#' }
#' All of these files must be present in the same directory for the constructor
#' \code{LDZipMatrix(prefix)} to succeed.
#'
#' @param prefix Character string. The file prefix of the LDZip matrix (without extensions).
#' @return An object of class \code{"LDZipMatrix"} (an external pointer to a C++ instance).
#'
#' @aliases LDZipMatrix LDZipMatrix-class
#'
#' @examples
#' \dontrun{
#' prefix <- file.path(system.file("extdata", package = "LDZipMatrix"), "g1k.chr22.ldzip")
#'
#' # Create a new LDZipMatrix handle
#' ld <- LDZipMatrix(prefix)
#'
#' # Query metadata
#' dim(ld)
#' stats(ld)
#'
#' # Print and summarize
#' print(ld)
#' summary(ld)
#' 
#' # Fetch a single entry (on-demand decompression)
#' val <- fetchLD(ld, 10, 10, type="PHASED_R")
#'
#' }
#'
#' @seealso \code{\link{LDZipMatrix}}, \code{\link{getLD}}, \code{\link{getNeighbors}}
#' 
#' @export
LDZipMatrix <- function(prefix) {
  ptr <- LDZipMatrix_rcpp(prefix)
  class(ptr) <- "LDZipMatrix"
  ptr
}
