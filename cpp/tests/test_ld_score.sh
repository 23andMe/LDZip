#!/bin/bash
set -e
OUT_DIR=output
rm -rf $OUT_DIR/*
source ./check_plink.sh

echo -e "\n\033[1;33m ----> Test LD Score \033[0m"

WINDOW_KB=100
PLINK_R2_THRESHOLD=0.01

echo "➤ Generating LD matrix with PLINK2..."
${PLINK2} \
  --pfile ../../assets/g1k \
  --ld-window-kb $WINDOW_KB \
  --ld-window-r2 $PLINK_R2_THRESHOLD \
  --r2-unphased ref-based cols=id,ref,alt \
  --out $OUT_DIR/plink_tabular > /dev/null

echo "➤ Compressing with ldzip..."
../bin/ldzip compress plinkTabular \
  --ld_file $OUT_DIR/plink_tabular.vcor \
  --snp_file ../../assets/g1k.pvar \
  --output_prefix $OUT_DIR/compressed \
  --bits 99 \
  --min 0 \
  --min_col UNPHASED_R2 > /dev/null

echo "➤ Running LD score calculation..."
../bin/ldzip ld-score \
  --input_prefix $OUT_DIR/compressed \
  --output $OUT_DIR/ldscores.txt \
  --window $WINDOW_KB \
  --stat UNPHASED_R2 > /dev/null

echo "➤ Validating LD scores with R..."
Rscript ../../scripts/check_ldscores.R \
  $OUT_DIR/plink_tabular.vcor \
  ../../assets/g1k.pvar \
  $OUT_DIR/ldscores.txt \
  $WINDOW_KB

echo "✅ LD score test passed!"
