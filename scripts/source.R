detect_ld_type <- function(path) {
  if (length(path) > 1) {
    if (all(grepl("\\.vcor$", path))) {
      return("vcor_list")
    } else if (all(grepl("\\.bin$", path))) {
      return("bin_list")
    } else {
      stop("Unknown list input type")
    }
  } else if (file.exists(paste0(path, "_i.bin"))) {
    return("compressed")
  } else if (grepl("\\.vcor$", path)) {
    return("vcor")
  } else if (grepl("\\.bin$", path)) {
    return("plink_bin")
  } else {
    stop("Unknown input type")
  }
}

read_ld_data <- function(path) {
  type <- detect_ld_type(path)
  if (type == "vcor") {
    df <- read.table(path, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
    stopifnot(ncol(df) == 7)
    return(as.numeric(df[[7]]))
  } else if (type == "plink_bin") {
    return(as.vector(read_plink_ld_matrix(path)))
  } else if (type == "compressed") {
    return(as.vector(read_dsC_encoded_full(path)))
  } else if (type == "bin_list") {
    mats <- lapply(path, read_plink_ld_matrix)
    return(as.vector(Matrix::bdiag(mats)))
  } else if (type == "vcor_list") {
    mats <- lapply(path, function(p) {
      df <- read.table(p, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
      stopifnot(ncol(df) == 7)
      as.numeric(df[[7]])
    })
  return(unlist(mats, use.names = FALSE))
  } else {
    stop("Unknown input type")
  }
}

compare_ld_headers <- function(orig_path, new_path) {

  orig_df <- read.table(orig_path, header = FALSE, sep = "\t", stringsAsFactors = FALSE, colClasses = "character")
  new_df  <- read.table(new_path,  header = FALSE, sep = "\t", stringsAsFactors = FALSE, colClasses = "character")

  stopifnot(ncol(orig_df) >= 6, ncol(new_df) >= 6)

  orig_keys <- orig_df[, 1:6]
  new_keys  <- new_df[, 1:6]

  if (!identical(orig_keys, new_keys)) {
    stop("First 6 columns of the input files do not match")
  }
}

remove_empty_ld_columns <- function(R) {
  i = 1
  while (i <= ncol(R) && sum(is.na(R[, i])) == nrow(R)) {i = i + 1}
  bad.ld <- is.na(R[, i])
  cat(Sys.time(), ": found [", sum(bad.ld), "] bad variants in R2 matrix")
  cat(Sys.time(), ": removing [", sum(bad.ld), "] bad variants from R2 matrix")
  R[!bad.ld, !bad.ld]
}

read_plink_ld_matrix = function(path){

  ld_file = paste0(path)
  snp_file = paste0(path, ".vars")
  N <- length(read.table(snp_file)[,1])
  R = matrix(readBin(paste0(ld_file), what="numeric", size=4, n=N**2), ncol=N)
  remove_empty_ld_columns(R)
}

write_filtered_ld_matrix = function(in_path, keep_file, out_path) {
  R <- read_plink_ld_matrix(in_path)

  keep_idx <- scan(keep_file, what = integer())
  keep_idx <- keep_idx + 1  # convert 0-based to 1-based

  Rf <- R[keep_idx, keep_idx, drop = FALSE]

  writeBin(as.numeric(Rf), out_path, size = 4)
}

read_tabular_ld_matrix <- function(path) {
  snp_file <- paste0(path, ".vars")
  ds  <- read.table(snp_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "")
  ids <- ds$ID
  N   <- length(ids)

  df  <- read.table(paste0(path, ".vcor"), header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "", check.names = FALSE)

  R <- matrix(0.0, nrow = N, ncol = N, dimnames = list(ids, ids))

  for (i in seq_len(nrow(df))) {
    a <- df$`#ID_A`[i]
    b <- df$ID_B[i]
    val <- df$UNPHASED_R[i]
    R[a, b] <- val
    R[b, a] <- val
  }

  diag(R) <- 1
  R
}

compare_values <- function(old_val, new_val) {
  # Metrics
  num_non_zero_old <- sum(old_val != 0.0)
  num_non_zero_new <- sum(new_val != 0.0)
  cor_full         <- cor(old_val, new_val)
  cor_nonzero      <- cor(old_val[old_val != 0], new_val[old_val != 0])
  max_diff         <- max(abs(old_val - new_val))

  list(
    num_non_zero_old = num_non_zero_old,
    num_non_zero_new = num_non_zero_new,
    cor_full         = cor_full,
    cor_nonzero      = cor_nonzero,
    max_diff         = max_diff
  )
}

get_threshold <- function(bits, min) {
  thresh <- read.table("../../../cpp/tests/assets/thresholds.tsv",
                       header = FALSE,
                       sep = "\t",
                       stringsAsFactors = FALSE)
  # assume file has columns: bits, min, threshold
  colnames(thresh) <- c("bits", "min", "threshold")

  val <- as.numeric(thresh$threshold[thresh$bits == bits & thresh$min == min])

  if (length(val) == 0) {
    stop(sprintf("No threshold found for bits=%s, min=%s", bits, min))
  }
  val
}
