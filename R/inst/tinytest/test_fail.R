library(LDZipMatrix)
library(tinytest)

prefix <- file.path(system.file("extdata", package = "LDZipMatrix"), "g1k.chr22.ldzip")
ld <- LDZipMatrix(prefix)
buildIndex(ld)

## -----------------------------
## empty inputs
## -----------------------------

expect_error(
  fetchLD(ld, integer(0), 1),
  "non-empty",
  info = "row empty"
)

expect_error(
  fetchLD(ld, 1, integer(0)),
  "non-empty",
  info = "col empty"
)

## -----------------------------
## duplicate entries
## -----------------------------

expect_error(
  fetchLD(ld, c(1, 1), 2),
  "duplicate",
  info = "duplicate row indices"
)

expect_error(
  fetchLD(ld, 1, c(2, 2)),
  "duplicate",
  info = "duplicate col indices"
)

expect_error(
  fetchLD(ld, c("rs587755077", "rs587755077"), c("rs587755077")),
  "duplicate",
  info = "duplicate row rsids"
)


## -----------------------------
## invalid numeric indices
## -----------------------------

expect_error(
  fetchLD(ld, 0, 1),
  "positive",
  info = "row index zero"
)

expect_error(
  fetchLD(ld, 1.5, 2),
  "integer",
  info = "row index non-integer"
)

## -----------------------------
## mixed row / col types
## -----------------------------

expect_error(
  fetchLD(ld, 1, "rs587755077"),
  "same type",
  info = "row numeric, col character"
)

expect_error(
  fetchLD(ld, "rs587755077", 1),
  "same type",
  info = "row character, col numeric"
)

## -----------------------------
## missing rsids
## -----------------------------
expect_error(
  fetchLD(ld, c("rs587755077", "rsFAKE123"), "rs587631919", type="PHASED_R"),
  "Missing",
  info = "missing rsid in row"
)
