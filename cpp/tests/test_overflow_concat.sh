#!/bin/bash
set -e
OUT_DIR=output
rm -rf $OUT_DIR/*
N=1000000
M=100
source ./check_plink.sh

echo -e "\n\033[1;33m ----> Test Overflow Tabular (Lossless) \033[0m"

echo "➤ Generating random matrix ..."
for chr in 20 21 22; do
  Rscript ../../scripts/generate_large_random_tabular.R ${N} output/test_tabular.${chr} ${M} > /dev/null
done

echo "➤ Compressing..."
for chr in 20 21 22; do
  ../bin/ldzip compress plinkTabular \
    --ld_file $OUT_DIR/test_tabular.${chr}.vcor \
    --snp_file $OUT_DIR/test_tabular.${chr}.vars \
    --output_prefix $OUT_DIR/compressed.${chr} \
    --min 0.0 \
    --bits 99 \
    --min_col UNPHASED_R > /dev/null
done


echo "➤ Concating..."
../bin/ldzip concat \
  --inputs $OUT_DIR/compressed.20 $OUT_DIR/compressed.21 $OUT_DIR/compressed.22 \
  --output_prefix $OUT_DIR/concat > /dev/null


echo "➤ Decompressing..."
../bin/ldzip decompress \
  --input_prefix $OUT_DIR/concat \
  --output_prefix $OUT_DIR/test_out \
  --type tabular


echo "➤ Comparing roundtrip..."
Rscript ../../scripts/corr.R \
    -orig $OUT_DIR/test_tabular.20.vcor $OUT_DIR/test_tabular.21.vcor $OUT_DIR/test_tabular.22.vcor\
    -new $OUT_DIR/test_out.vcor > /dev/null
./test_identical.sh $OUT_DIR/test_out_stats.json
