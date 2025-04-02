#!/bin/bash
set -e
OUT_DIR=output
rm -rf $OUT_DIR/*
N=1000
source ./check_plink.sh

echo -e "\n\033[1;33m ----> Test Filtering PLINK Tabular Compression (Lossless) \033[0m"


echo "➤ Generating random plink matrix ..."
${PLINK2} \
  --pfile ../../assets/g1k \
  --geno 0 \
  --write-snplist \
  --out output/temp_filtered > /dev/null
sort -R output/temp_filtered.snplist | head -n "$N" > output/top.snps
sort -R output/top.snps | head -n $((N/2)) > output/filtered.snps
grep -F -n -f output/filtered.snps ../../assets/g1k.pvar | cut -d: -f1 | awk '{print $1 - 2}' | sort -n > output/filtered.indices

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

echo "➤ Filtering compressed ..."
../bin/ldzip filter \
  --input_prefix $OUT_DIR/compressed \
  --keep $OUT_DIR/filtered.indices \
  --output_prefix $OUT_DIR/compressed_filtered  > /dev/null

echo "➤ Decompressing..."
../bin/ldzip decompress \
  --input_prefix $OUT_DIR/compressed_filtered \
  --output_prefix $OUT_DIR/test_out \
  --type tabular > /dev/null

echo "➤ Creating uncompressed filtered using plink..."
${PLINK2} \
  --pfile ../../assets/g1k \
  --extract  output/filtered.snps \
  --ld-window-kb 100 \
  --ld-window-r2 0.000001 \
  -r-unphased ref-based cols=id,ref,alt \
  --out $OUT_DIR/test_in_filtered > /dev/null > /dev/null

echo "➤ Comparing roundtrip..."
Rscript ../../scripts/corr.R \
    -orig $OUT_DIR/test_in_filtered.vcor \
    -new $OUT_DIR/test_out.vcor > /dev/null

./test_identical.sh $OUT_DIR/test_out_stats.json
