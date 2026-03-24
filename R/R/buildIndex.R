#' Build an index for querying by variant ID from an LDZipMatrix object
#'
#' This function builds an index on variant identifiers from the 
#' associated \code{.vars.txt} file of the \code{LDZipMatrix} object. 
#' Optionally one can pass their own variant file. 
#' 
#' The index allows efficient mapping from variant IDs to internal 
#' integer indices, which can then be used with \code{fetchLD()} 
#' or \code{getNeighbors()}.
#'
#' @param ld An external pointer to an \code{LDZipMatrix} object,
#'   typically created with \code{LDZipMatrix()}.
#' @param variant_file optional path to the variant annotation file (e.g. a 
#'   \code{.pvar} file) containing the variant information 
#'   (\code{CHROM}, \code{POS}, \code{ID}, \code{REF}, \code{ALT}). A standard PLINK2 \code{.pvar} file works here.
#' @return
#'   \code{NULL}. The index is written to the SQLite database (\code{$prefix.sqlite})
#' @examples
#' \dontrun{
#' # Open an LDZip matrix
#' prefix <- file.path(system.file("extdata", package = "LDZipMatrix"), "g1k.chr22.ldzip")
#' ld <- LDZipMatrix(prefix)
#' 
#' # Build an index
#' idx <- buildIndex(ld)
#'
#' # Fetch LD values by variant IDs
#' fetchLD(ld, "rs587755077", c("rs587631919", "rs587661542"), type = "PHASED_R")
#' }
#'
#' @export
buildIndex <- function(ld, variant_file) {

    if (missing(variant_file))
        variant_file <- paste(LDZipMatrix_get_prefix_rcpp(ld), "vars", "txt", sep = ".")
    if (!file.exists(variant_file)) {
        stop(sprintf("Variant file not found: %s", variant_file))
    }
    db_file <- paste(LDZipMatrix_get_prefix_rcpp(ld), "sqlite", sep = ".")

    message(date(), ": reading variant file")
    variants <- read.table(
        variant_file,
        comment.char = "#",
        header = FALSE,
        stringsAsFactors = FALSE
    )[, 1:5]
    colnames(variants) <- c("CHROM", "POS", "ID", "REF", "ALT")

    temp <- parse_variant_identifier(variants$ID)
    variants$RSLEV <- temp$fmt
    variants$RSNUM <- as.character(temp$num)   # safe for insertion
    variants$ID <- NULL

    con <- DBI::dbConnect(RSQLite::SQLite(), db_file)
    on.exit(DBI::dbDisconnect(con), add = TRUE)

    message(date(), ": creating schema")
    DBI::dbExecute(con, "DROP TABLE IF EXISTS variants")
    DBI::dbExecute(con, "
        CREATE TABLE variants (
          CHROM  INTEGER,
          POS    INTEGER,
          RSNUM  INTEGER,
          RSLEV  TEXT,
          REF    TEXT,
          ALT    TEXT
        )
    ")

    message(date(), ": writing variant table to database")
    DBI::dbExecute(con, "BEGIN TRANSACTION")
    DBI::dbWriteTable(con, "variants", variants, append = TRUE)
    DBI::dbExecute(con, "COMMIT")

    message(date(), ": building UNIQUE composite index on RSNUM, RSLEV")
    DBI::dbExecute(con, "
        CREATE UNIQUE INDEX IF NOT EXISTS idx_variants_rsnum_rslev
        ON variants(RSNUM, RSLEV)
    ")

    message(date(), ": building composite index on genomic region (CHROM, POS)")
    DBI::dbExecute(con, "
        CREATE INDEX IF NOT EXISTS idx_variants_region
        ON variants(CHROM, POS)
    ")

    invisible(db_file)
}
