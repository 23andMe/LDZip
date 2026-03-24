.check_fetchLD_inputs <- function(row, col, pairwise) {

  if (!length(row) || !length(col))
    stop("`row` and `col` must be non-empty")

	if (is.numeric(row) != is.numeric(col))
  	stop("`row` and `col` must be of the same type (both numeric or both character)")

  if (pairwise) {
    if (length(row) != length(col))
      stop("`pairwise=TRUE` requires `row` and `col` to have equal length")
  } else {
    if (anyDuplicated(row))
      stop("`row` contains duplicate entries")
    if (anyDuplicated(col))
      stop("`col` contains duplicate entries")
  }

  if (is.numeric(row) && any(row < 1 | row %% 1 != 0))
    stop("Numeric `row` indices must be positive integers")
  if (is.numeric(col) && any(col < 1 | col %% 1 != 0))
    stop("Numeric `col` indices must be positive integers")

  invisible(TRUE)
}

#' Fetch LD values from an LDZipMatrix object
#'
#' This function retrieves linkage disequilibrium (LD) values from a
#' compressed \code{LDZipMatrix} object. Queries can be made either by
#' integer indices (1-based) or by variant identifiers (IDs). If both
#' indices and variant IDs are provided, variant IDs take priority.
#'
#' @param ld An external pointer to an \code{LDZipMatrix} object,
#'   typically created with \code{LDZipMatrix()}.
#' @param row,col An Integer (1-based) or character scalar/vector of variant indices/IDs
#'   specifying the row(s)/column(s) to query.
#' @param types a character scalar/vector of types of linkage statistics to 
#'   return. Choices are \code{"PHASED_R"}, \code{"UNPHASED_R"}, \code{"PHASED_R2"}, \code{"UNPHASED_R2"}, \code{"D"}, \code{"DPRIME"}.
#' @param pairwise Logical flag controlling the return type:
#'   \itemize{
#'     \item \code{TRUE}: return pairwise LD values, where each element of
#'       \code{row} is matched with the corresponding element of \code{col}
#'       (both must have equal length).
#'     \item \code{FALSE}: return the full LD submatrix spanning all
#'       combinations of \code{row} and \code{col}.
#'   }
#' @param simplify Logical flag. If \code{TRUE}, the return value is simplified
#'   according to the shape of the query and the number of statistics:
#'   \itemize{
#'     \item Single statistic, (\code{1 x 1}), (\code{1 x N}) or (\code{N x 1}): a \code{vector}.
#'     \item Single statistic, (\code{N x N}): a \code{matrix}.
#'     \item Multiple statistics, (\code{1 x 1}), (\code{1 x N}) or (\code{N x 1}): a \code{data.frame}.
#'     \item Multiple statistics, (\code{N x N}): a \code{list} of \code{matrix}.
#'   }
#'   If \code{FALSE}, results are always returned as a list of matrices, one per statistic.
#'   Note: This parameter is ignored when \code{pairwise=TRUE}.
#'
#' @return
#'   When \code{pairwise=FALSE} (default):
#'   \itemize{
#'     \item A numeric vector if either \code{row} or \code{col} is a vector
#'       (the other scalar), and only one statistic is requested.
#'     \item A numeric matrix if both \code{row} and \code{col} are vectors and
#'       a single statistic is requested.
#'     \item A \code{data.frame} if either \code{row} or \code{col} is a vector
#'       (the other scalar), and multiple statistics are requested.
#'     \item A list of matrices if both \code{row} and \code{col} are vectors
#'       and multiple statistics are requested (or if \code{simplify=FALSE})
#'   }
#'   When \code{pairwise=TRUE}:
#'   \itemize{
#'     \item Always returns a \code{data.frame} with one column per requested statistic type,
#'       one row per pair, and additional columns \code{var1} and \code{var2} to identify the pairs.
#'       Results are returned in the same order as the input. Duplicate pairs are allowed.
#'   }
#'
#' @note
#' \itemize{
#'   \item In order to use variant identifiers (instead of indices), please index the \code{LDZipMatrix} object using \code{buildIndex(ld)}
#'   \item When \code{pairwise=FALSE}: Variants in the returned object are sorted in order of their variant indices, and not the same order as the input in \code{row} or \code{col}
#'   \item When \code{pairwise=TRUE}: Results are returned in the same order as the input. Duplicate pairs are allowed.
#' }
#' @examples
#' \dontrun{
#' # Open an LDZip matrix
#' prefix <- file.path(system.file("extdata", package = "LDZipMatrix"), "g1k.chr22.ldzip")
#' ld <- LDZipMatrix(prefix)
#' 
#' # 1 x 1 query, single statistic → scalar
#' fetchLD(ld, 2, 191, types = "PHASED_R")
#'
#'  # 1 x 1 query, multiple statistics → data.frame
#' fetchLD(ld, 2, 191, types = c("PHASED_R", "DPRIME"))
#'
#' # 1 x N query, single statistic → vector
#' fetchLD(ld, 2, c(2,129,191,307), types = "PHASED_R")
#' 
#' # 1 x N query, multiple statistics → data.frame
#' fetchLD(ld, 2, c(2,129,191,307), types =  c("PHASED_R", "DPRIME"))
#'
#' # N x N query, single statistic → matrix
#' fetchLD(ld, 1:5, 1:5, types = "PHASED_R")
#'
#' # N x N query, multiple statistics → list of matrices
#' fetchLD(ld, 1:5, 1:5, types = c("PHASED_R", "DPRIME"))
#' 
#' ## Query with variant identifiers
#' # Build index (if it doesn't exist already)
#' buildIndex(ld)
#' 
#' # 1 x 1 query, single statistic → scalar
#' fetchLD(ld, "rs587755077", "rs587631919", types = "PHASED_R")
#'
#' # Define a small set of variant IDs
#' vars <- c("rs587725733", "rs587631919", "rs587661542", "rs587718290", "rs587729333")
#' 
#' # 1 x N query, single statistic → vector
#' fetchLD(ld, "rs587755077", vars, types = "PHASED_R")
#'
#' # 1 x N query, multiple statistics → data.frame
#' fetchLD(ld, "rs587755077", vars, types = c("PHASED_R", "DPRIME"))
#'
#' # N x N query, single statistic → matrix
#' fetchLD(ld, vars, vars, types = "PHASED_R")
#'
#' # N x N query, multiple statistics → list of matrices
#' fetchLD(ld, vars, vars, types = c("PHASED_R", "DPRIME"))
#'
#' ## Pairwise queries (always returns data.frame)
#' # Fetch LD for specific pairs (row[i] with col[i])
#' rows <- c(1, 5, 10)
#' cols <- c(2, 8, 15)
#' fetchLD(ld, rows, cols, types = "PHASED_R", pairwise = TRUE)
#'
#' # Pairwise with variant IDs
#' row_vars <- c("rs587725733", "rs587631919", "rs587661542")
#' col_vars <- c("rs587755077", "rs587718290", "rs587729333")
#' fetchLD(ld, row_vars, col_vars, types = "PHASED_R", pairwise = TRUE)
#'
#' # Pairwise with multiple statistics (returns data.frame with multiple columns)
#' fetchLD(ld, rows, cols, types = c("PHASED_R", "DPRIME"), pairwise = TRUE)
#' }
#'
#' @export
fetchLD <- function(ld, row, col, types=c("UNPHASED_R"), pairwise=FALSE, simplify=TRUE) {
  
	.check_fetchLD_inputs(row, col, pairwise)
	rowNames = row
	colNames = col

	if (is.character(row) || is.character(col)) {
		variant_db_file <- paste(LDZipMatrix_get_prefix_rcpp(ld), "sqlite", sep=".")
		if (is.character(row)) { r <- get_indices_by_rsid(variant_db_file, rsids=row); row_idx <- r$idx; rowNames <- r$rsid } else { row_idx <- row }
		if (is.character(col)) { c <- get_indices_by_rsid(variant_db_file, rsids=col); col_idx <- c$idx; colNames <- c$rsid } else { col_idx <- col }
	} else {
		row_idx <- row
		col_idx <- col
	}

	# Sort or validate based on mode
	if (!pairwise) {
		o<-order(row_idx); row_idx<-row_idx[o]; rowNames<-rowNames[o]
		o<-order(col_idx); col_idx<-col_idx[o]; colNames<-colNames[o]
	}

	# Handle pairwise case
	if (pairwise) {
		# Get pairwise values for each stat type (for each col, the rows should be sorted)
		vals <- lapply(types, function(type) LDZipMatrix:::LDZipMatrix_get_pairswise_rcpp(ld, as.integer(col_idx - 1L), as.integer(row_idx - 1L), type))

		# Build data.frame with types as columns
		res <- as.data.frame(vals, col.names = types)
		res$var1 <- rowNames
		res$var2 <- colNames
		rownames(res) <- NULL

		return(res)
	}

	is_range <- function(x) length(x) <= 1L || all(diff(x) == 1L)

	# One by One
	if (length(row_idx) == 1L && length(col_idx) == 1L) {

		vals <- lapply(types, function(type) LDZipMatrix_getValue_rcpp(ld, row_idx - 1L, col_idx - 1L, type))
		if (simplify) 
			res <- if (length(types) == 1L) vals[[1]] else as.data.frame(vals, col.names = types)
		else 
			res <- setNames(lapply(vals, function(v) matrix(v, length(row_idx), length(col_idx), dimnames = list(rowNames, colNames))), types)
		
	} else if (length(row_idx) == 1L) {

		vals <- lapply(types, function(type) {
			if (is_range(col_idx)) 
				LDZipMatrix:::LDZipMatrix_get_column_range_rcpp(ld, row_idx - 1L, col_idx[1L] - 1L, col_idx[length(col_idx)] - 1L, type)
			else 
				LDZipMatrix:::LDZipMatrix_get_column_indices_rcpp(ld, row_idx - 1L, col_idx - 1L, type)
		})
		if (simplify) {
			if (length(types) == 1L) {
				res <- vals[[1L]]
				names(res) <- colNames
			} else {
			res <- as.data.frame(vals, col.names = types, row.names = colNames)
			}
		} else {
			res <- setNames(lapply(vals, function(v) matrix(v, length(row_idx), length(col_idx), dimnames = list(rowNames, colNames))), types)
		}

	} else if (length(col_idx) == 1L) {

		vals <- lapply(types, function(type) {
			if (is_range(row_idx)) 
				LDZipMatrix:::LDZipMatrix_get_column_range_rcpp(ld, col_idx - 1L, row_idx[1L] - 1L, row_idx[length(row_idx)] - 1L, type)
			else 
				LDZipMatrix:::LDZipMatrix_get_column_indices_rcpp(ld, col_idx - 1L, row_idx - 1L, type)
		})
		if (simplify) {
			if (length(types) == 1L) {
				res <- vals[[1L]]
				names(res) <- rowNames
			} else {
				res <- as.data.frame(vals, col.names = types, row.names = rowNames)
			}
		} else {
			res <- setNames(lapply(vals, function(v) matrix(v, length(row_idx), length(col_idx), dimnames = list(rowNames, colNames))), types)
		}

	} else {

		vals <- lapply(types, function(type) {
			if (is_range(row_idx) && is_range(col_idx)) 
				LDZipMatrix:::LDZipMatrix_get_submatrix_range_rcpp(ld, col_idx[1L] - 1L, col_idx[length(col_idx)] - 1L,row_idx[1L] - 1L, row_idx[length(row_idx)] - 1L, type)
			else
				LDZipMatrix:::LDZipMatrix_get_submatrix_indices_rcpp(ld, col_idx - 1L, row_idx - 1L, type)
		})
		res <- setNames(lapply(vals, function(m) { dimnames(m) <- list(rowNames, colNames); m }), types)
		if (simplify && length(types) == 1L) {
			res <- res[[1L]]
		}
	}
	res

	}
