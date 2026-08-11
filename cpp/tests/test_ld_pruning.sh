#!/bin/bash
set -e
OUT_DIR=output
rm -rf $OUT_DIR/*
source ./check_plink.sh

echo -e "\n\033[1;33m ----> Test LD Pruning \033[0m"

WINDOW_KB=100
R_THRESHOLD=0.7
MIN=0.1

echo "➤ Generating LD matrix with PLINK2..."
min_r2=$(echo "$MIN * $MIN" | bc -l)
${PLINK2} \
  --pfile ../../assets/g1k \
  --ld-window-kb $WINDOW_KB \
  --ld-window-r2 $min_r2 \
  --r-unphased ref-based cols=id,ref,alt \
  --out $OUT_DIR/plink_tabular > /dev/null

echo "➤ Compressing with ldzip..."
../bin/ldzip compress plinkTabular \
  --ld_file $OUT_DIR/plink_tabular.vcor \
  --snp_file ../../assets/g1k.pvar \
  --output_prefix $OUT_DIR/compressed \
  --bits 99 \
  --min $MIN \
  --min_col UNPHASED_R > /dev/null

echo "➤ Running LD pruning..."
../bin/ldzip prune \
  --input_prefix $OUT_DIR/compressed \
  --output_prefix $OUT_DIR/pruned \
  --window $WINDOW_KB \
  --threshold $R_THRESHOLD \
  --stat UNPHASED_R > /dev/null

echo "➤ Validating output..."

TOTAL_VARIANTS=$(tail -n +2 ../../assets/g1k.pvar | wc -l)
KEPT_COUNT=$(wc -l < $OUT_DIR/pruned.prune.in)
REMOVED_COUNT=$(wc -l < $OUT_DIR/pruned.prune.out)

echo "   - Total variants: $TOTAL_VARIANTS"
echo "   - Kept variants:  $KEPT_COUNT"
echo "   - Removed variants: $REMOVED_COUNT"

# Check that kept + removed = total
SUM=$((KEPT_COUNT + REMOVED_COUNT))
if [ $SUM -ne $TOTAL_VARIANTS ]; then
    echo "❌ Kept ($KEPT_COUNT) + Removed ($REMOVED_COUNT) != Total ($TOTAL_VARIANTS)"
    exit 1
fi

# Check no duplicates
UNIQUE_KEPT=$(sort -u $OUT_DIR/pruned.prune.in | wc -l)
if [ $UNIQUE_KEPT -ne $KEPT_COUNT ]; then
    echo "❌ Duplicate variants in .prune.in"
    exit 1
fi

UNIQUE_REMOVED=$(sort -u $OUT_DIR/pruned.prune.out | wc -l)
if [ $UNIQUE_REMOVED -ne $REMOVED_COUNT ]; then
    echo "❌ Duplicate variants in .prune.out"
    exit 1
fi

# Check no overlap
OVERLAP=$(cat $OUT_DIR/pruned.prune.in $OUT_DIR/pruned.prune.out | sort | uniq -d | wc -l)
if [ $OVERLAP -gt 0 ]; then
    echo "❌ Overlap between .prune.in and .prune.out"
    exit 1
fi

# Check that some variants were pruned
if [ $REMOVED_COUNT -eq 0 ]; then
    echo "❌ No variants were pruned"
    exit 1
fi

echo "✅ All checks passed!"
echo "   - Pruning rate: $(echo "scale=2; 100 * $REMOVED_COUNT / $TOTAL_VARIANTS" | bc)%"

echo "➤ Validating LD thresholds with R..."
Rscript ../../scripts/check_pruned.R \
  $OUT_DIR/plink_tabular.vcor \
  ../../assets/g1k.pvar \
  $OUT_DIR/pruned.prune.in \
  $OUT_DIR/pruned.prune.out \
  $R_THRESHOLD \
  $WINDOW_KB
