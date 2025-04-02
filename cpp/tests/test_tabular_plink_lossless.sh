#!/bin/bash
set -e
OUT_DIR=output
rm -rf $OUT_DIR/*
N=1000
source ./check_plink.sh

echo -e "\n\033[1;33m ----> Test PLINK Tabular Compression (Lossless) \033[0m"


echo "➤ Generating random plink matrix ..."
${PLINK2} \
  --pfile ../../assets/g1k \
  --geno 0 \
  --write-snplist \
  --out output/temp_filtered > /dev/null
sort -R output/temp_filtered.snplist | head -n "$N" > output/top.snps
${PLINK2} \
  --pfile ../../assets/g1k \
  --extract  output/top.snps \
  --ld-window-kb 100 \
  --ld-window-r2 0.000001 \
  -r-unphased ref-based cols=id,ref,alt \
  --out $OUT_DIR/plink_tabular > /dev/null


echo "➤ Compressing..."
../bin/ldzip compress plinkTabular \
  --ld_file $OUT_DIR/plink_tabular.vcor \
  --snp_file ../../assets/g1k.pvar \
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
    -orig $OUT_DIR/plink_tabular.vcor \
    -new $OUT_DIR/test_out.vcor > /dev/null

./test_identical.sh $OUT_DIR/test_out_stats.json
