library(LDZipMatrix)
library(tinytest)

prefix <- file.path(system.file("extdata", package = "LDZipMatrix"), "g1k.chr22.ldzip")
ld <- LDZipMatrix(prefix)
buildIndex(ld)

# define test sets
idx1  <- 2
idxN  <- c(5, 1, 3, 4, 2)
rsid1 <- "rs587755077"
rsidN <- c("rs587718290",
           "rs587725733",
           "rs587661542",
           "rs587729333",
           "rs587631919")

# actual order hard coded
canon_rsid <- c("rs587725733",
                "rs587631919",
                "rs587661542",
                "rs587718290",
                "rs587729333")

for (simplify in c(TRUE, FALSE)) {
  for (use_rsid in c(FALSE, TRUE)) {
    for (multi in c(FALSE, TRUE)) {

      one   <- if (use_rsid) rsid1 else idx1
      many  <- if (use_rsid) rsidN else idxN
      types <- if (multi) c("PHASED_R", "DPRIME") else "PHASED_R"

      msg <- paste("order check:",
                   "simplify =", simplify,
                   "use_rsid =", use_rsid,
                   "multi =", multi)

      out <- fetchLD(ld, one, many, type = types, simplify = simplify)

      # extract returned names
      if (simplify && !multi) {
        returned <- names(out)
      } else if (simplify && multi) {
        returned <- rownames(out)
      } else {
        returned <- colnames(out[[1]])
      }

      # expected canonical order
      expected <- if (use_rsid) canon_rsid else as.character(sort(idxN))

      expect_equal(returned, expected, info = msg)
    }
  }
}
