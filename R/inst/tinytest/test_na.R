library(LDZipMatrix)
library(tinytest)

ld <- LDZipMatrix("../../tests/data/compressed_na")

set.seed(1)

test_cases <- replicate(
  10,
  {
    i <- sample(1:1000, 1)
    j <- sample(i:1000, 1)
    i:j
  },
  simplify = FALSE
)

for (tc in test_cases) {
  mat <- fetchLD(ld, tc, tc, type="D")

  expect_true(
    all(is.na(diag(mat))),
    info = paste0("Diagonal not NA for range [", min(tc), ", ", max(tc), "]")
  )
}
