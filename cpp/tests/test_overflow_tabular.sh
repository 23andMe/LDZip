#!/bin/bash
set -e
OUT_DIR=output
rm -rf $OUT_DIR/*
N=1000
source ./check_plink.sh

echo -e "\n\033[1;33m ----> Test Overflow Tabular (Lossless) \033[0m"


echo "➤ Generating random matrix ..."
Rscript ../../scripts/generate_large_random_tabular.R 100000 output/test_tabular 20 > /dev/null


echo "➤ Compressing..."
../bin/ldzip compress plinkTabular \
  --ld_file $OUT_DIR/test_tabular.vcor \
  --snp_file $OUT_DIR/test_tabular.vars \
  --output_prefix $OUT_DIR/compressed \
  --min 0.0 \
  --bits 99 \
  --min_col UNPHASED_R > /dev/null


echo "➤ Decompressing..."
../bin/ldzip decompress \
  --input_prefix $OUT_DIR/compressed \
  --output_prefix $OUT_DIR/test_out \
  --type tabular

echo "➤ Comparing roundtrip..."
Rscript ../../scripts/corr.R \
    -orig $OUT_DIR/test_tabular.vcor \
    -new $OUT_DIR/test_out.vcor > /dev/null
./test_identical.sh $OUT_DIR/test_out_stats.json
