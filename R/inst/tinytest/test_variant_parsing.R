library(LDZipMatrix)
library(tinytest)


check_round_trip <- function(ids) {
  parsed <- LDZipMatrix:::parse_variant_identifier(ids)
  rejoined <- LDZipMatrix:::join_variant_identifier(parsed$num, parsed$fmt)
  all(ids == rejoined)
}


rsids = c(
    "rs123",
    "rs123_A_T",
    "rs123:A:T",
    "rs123456_chr17_A_T",
    "rs123456789993239_chr17_A_T",
    "chr17_100010000_rs123456789993239_A_T",
    "rs42",
    "chr5_98765_A_G",
    "rs1404283754",
    "rs1307108837:G:A",
    "rs1384315640",
    "rs1486704209",
    "rs1231210093:C:T",
    "rs1290631452",
    "rs1476192063",
    "rs1323278608",
    "chr1:52471:C:A",
    "rs1167114720",
    "chr1:48893:C:G",
    "rs1161466707",
    "chr1:60855:T:G",
    "rs192890528"
)

expect_true(check_round_trip(rsids))
expect_error(LDZipMatrix:::parse_variant_identifier("rs"))
expect_error(LDZipMatrix:::parse_variant_identifier(c("rs", "rsABC")))
