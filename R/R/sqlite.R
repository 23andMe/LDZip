get_query_sql = function(rsids)
{
  temp <- parse_variant_identifier(rsids)
  query_df <- data.frame(
    RSNUM = as.character(temp$num),
    RSLEV = temp$fmt,
    stringsAsFactors = FALSE
  )

  # Build UNION ALL SELECT subquery, for some reason I couldn't get the usual
  # in (...) to work for row values/tuples
  selects <- apply(query_df, 1, function(row) {
    sprintf("SELECT %s AS RSNUM, '%s' AS RSLEV", row[["RSNUM"]], row[["RSLEV"]])
  })
  union_sql <- paste(selects, collapse = " UNION ALL ")
}

parse_variant_identifier <- function(x, re = "^(.*\\D)(\\d+)(\\D*)$") {
  has_num <- grepl("\\d+", x)  # pre-check
  if (!all(has_num)) {
    bad <- x[!has_num]
    stop(sprintf("No numeric part found in identifier(s): %s", paste(bad, collapse=", ")))
  }
  num <- bit64::as.integer64(sub(re, "\\2", x, perl=TRUE))
  fmt <- sub(re, "\\1%\\3", x, perl=TRUE)
  list(num=num, fmt=fmt)
}

join_variant_identifier <- function(num, fmt){
  mapply(function(n,f) sub("%", as.character(n), f, fixed=TRUE), num, fmt, USE.NAMES=FALSE)
}

get_indices_by_rsid <- function(db_file, rsids) {
  if (!file.exists(db_file)) stop(sprintf("Database file not found: %s", db_file))

  con <- DBI::dbConnect(RSQLite::SQLite(), db_file)
  on.exit(DBI::dbDisconnect(con), add = TRUE)

  union_sql = get_query_sql(rsids)
  sql <- sprintf("
    SELECT v.rowid AS idx, v.RSNUM, v.RSLEV
    FROM variants v
    JOIN ( %s ) AS bar
      ON v.RSNUM = bar.RSNUM AND v.RSLEV = bar.RSLEV
  ", union_sql)

  res <- DBI::dbGetQuery(con, sql)
  if (nrow(res) == 0) stop('No matching rsIDs found in variants table')

  # rebuild identifiers
  res$rsid <- join_variant_identifier(res$RSNUM, res$RSLEV)
  if (length(m <- setdiff(rsids, res$rsid))) stop("Missing rsIDs: ", paste(m, collapse = ", "))
  res[order(res$idx), c('idx','rsid')]
}

get_indices_by_region <- function(db_file, region) {
  if (!file.exists(db_file)) stop(sprintf("Database file not found: %s", db_file))
  stopifnot(all(c("chrom", "start", "end") %in% names(region)))

  con <- DBI::dbConnect(RSQLite::SQLite(), db_file)
  on.exit(DBI::dbDisconnect(con), add = TRUE)

  query <- sprintf(
    "SELECT rowid AS idx
     FROM variants
     WHERE CHROM = %s AND POS BETWEEN %d AND %d",
    region$chrom, region$start, region$end
  )
  res <- DBI::dbGetQuery(con, query)
  if (nrow(res) == 0) stop('No matching rsIDs found in specified region')

  sort(res$idx)
}

get_rsids_by_region <- function(db_file, region) {
  if (!file.exists(db_file)) stop(sprintf("Database file not found: %s", db_file))
  stopifnot(all(c("chrom", "start", "end") %in% names(region)))

  con <- DBI::dbConnect(RSQLite::SQLite(), db_file)
  on.exit(DBI::dbDisconnect(con), add = TRUE)

  sql <- sprintf("
    SELECT rowid AS idx, RSNUM, RSLEV
    FROM variants
    WHERE CHROM = %s AND POS BETWEEN %d AND %d
  ", region$chrom, region$start, region$end)

  res <- DBI::dbGetQuery(con, sql)
  if (nrow(res) == 0) stop('No matching rsIDs found in specified region')

  res$rsid <- join_variant_identifier(res$RSNUM, res$RSLEV)
  res <- res[order(res$idx), ]

  list(idx   = res$idx, rsids = res$rsid)
}

get_cpra_from_rsids <- function(db_file, rsids) {
  if (!file.exists(db_file)) stop(sprintf("Database file not found: %s", db_file))

  con <- DBI::dbConnect(RSQLite::SQLite(), db_file)
  on.exit(DBI::dbDisconnect(con), add = TRUE)

  union_sql = get_query_sql(rsids)
  sql <- sprintf("
    SELECT v.rowid AS idx, v.CHROM, v.POS, v.RSNUM, v.RSLEV
    FROM variants v
    JOIN ( %s ) AS q
      ON v.RSNUM = q.RSNUM AND v.RSLEV = q.RSLEV
  ", union_sql)

  res <- DBI::dbGetQuery(con, sql)
  if (nrow(res) == 0) stop('No matching rsIDs found in variants table')

  # rebuild identifiers
  res$ID <- join_variant_identifier(res$RSNUM, res$RSLEV)
  res[order(res$idx), c("idx", "CHROM", "POS", "ID")]
}

get_cpra_from_indices <- function(db_file, indices) {
  if (!file.exists(db_file)) stop(sprintf("Database file not found: %s", db_file))
  
  con <- DBI::dbConnect(RSQLite::SQLite(), db_file)
  on.exit(DBI::dbDisconnect(con), add = TRUE)

  sql <- sprintf(
    "SELECT rowid AS idx, CHROM, POS, RSNUM, RSLEV
     FROM variants 
     WHERE rowid IN (%s)",
    paste0(indices, collapse = ",")
  )
  res <- DBI::dbGetQuery(con, sql)

  # rebuild identifiers
  res$ID <- join_variant_identifier(res$RSNUM, res$RSLEV)
  res[order(res$idx), c("idx", "CHROM", "POS", "ID")]
}
