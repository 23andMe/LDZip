#!/bin/bash
set -e
OUT_DIR=output
rm -rf $OUT_DIR/*
N=1000
TYPE=PHASED_R
source ./check_plink.sh

echo -e "\n\033[1;33m ----> Test PLINK Compression (Lossless) \033[0m"


echo "➤ Generating random plink matrix ..."
${PLINK2} \
  --pfile ../../assets/g1k \
  --write-snplist \
  --out output/temp_filtered > /dev/null
sort -R output/temp_filtered.snplist | head -n "$N" > output/top.snps
${PLINK2} \
  --pfile ../../assets/g1k \
  --extract  output/top.snps \
  --snps-only \
  --r-unphased bin4 ref-based \
  --out output/plink_matrix > /dev/null


echo "➤ Compressing with ldzip using type [$TYPE]..."
../bin/ldzip compress plinkSquare \
  --ld_file $OUT_DIR/plink_matrix.unphased.vcor1.bin \
  --snp_file $OUT_DIR/plink_matrix.unphased.vcor1.bin.vars \
  --output_prefix $OUT_DIR/compressed \
  --min 0.0 \
  --bits 99 \
  --type $TYPE > /dev/null


echo "➤ Decompressing with ldzip..."
../bin/ldzip decompress \
  --input_prefix $OUT_DIR/compressed \
  --output_prefix $OUT_DIR/test_out


echo "➤ Comparing roundtrip..."
mv $OUT_DIR/test_out.$TYPE.bin $OUT_DIR/test_out.bin
Rscript ../../scripts/corr.R \
    -orig $OUT_DIR/plink_matrix.unphased.vcor1.bin \
    -new $OUT_DIR/test_out.bin > /dev/null


./test_identical.sh $OUT_DIR/test_out_stats.json

