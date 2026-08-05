#!/bin/bash
set -e
OUT_DIR=output
rm -rf $OUT_DIR/*
source ./check_plink.sh

echo -e "\n\033[1;33m ----> Test Find Tag Variants \033[0m"

THRESHOLD=0.8
MIN=0.01

echo "➤ Generating LD matrix with PLINK2..."
min_r2=$(echo "$MIN * $MIN" | bc -l)
${PLINK2} \
  --pfile ../../assets/g1k \
  --ld-window-kb 100 \
  --ld-window-r2 $min_r2 \
  --r-unphased ref-based cols=id,ref,alt \
  --out $OUT_DIR/plink_tabular > /dev/null

cp ../../assets/g1k.pvar $OUT_DIR/plink_tabular.vcor.vars.txt

echo "➤ Compressing with ldzip..."
../bin/ldzip compress plinkTabular \
  --ld_file $OUT_DIR/plink_tabular.vcor \
  --snp_file ../../assets/g1k.pvar \
  --output_prefix $OUT_DIR/compressed \
  --bits 99 \
  --min $MIN \
  --min_col UNPHASED_R > /dev/null

echo "➤ Creating test variant list (random 10 variants)..."
shuf -i 0-999 -n 10 | sort -n > $OUT_DIR/variants.txt

echo "➤ Running find-tag-variants..."
../bin/ldzip find-tag-variants \
  --input_prefix $OUT_DIR/compressed \
  --variants $OUT_DIR/variants.txt \
  --output $OUT_DIR/tags.txt \
  --threshold $THRESHOLD \
  --stat UNPHASED_R > /dev/null

echo "➤ Checking output file..."
if [ ! -s "$OUT_DIR/tags.txt" ]; then
    echo "❌ Missing or empty output file: $OUT_DIR/tags.txt"
    exit 1
fi

# Check that output has header
HEADER=$(head -1 $OUT_DIR/tags.txt)
if [ "$HEADER" != "variant	tag_variant" ]; then
    echo "❌ Invalid header in output file. Expected 'variant	tag_variant', got '$HEADER'"
    exit 1
fi

# Count lines (excluding header)
TAG_COUNT=$(tail -n +2 $OUT_DIR/tags.txt | wc -l)
echo "➤ Found $TAG_COUNT tag variants"

# Check that no variant tags itself
SELF_TAGS=$(tail -n +2 $OUT_DIR/tags.txt | awk '$1 == $2' | wc -l)
if [ $SELF_TAGS -gt 0 ]; then
    echo "❌ Found $SELF_TAGS self-tags (variants tagging themselves)"
    tail -n +2 $OUT_DIR/tags.txt | awk '$1 == $2' | head -5
    exit 1
fi

# Check that all variant indices are in input list
INPUT_VARIANTS=$(cat $OUT_DIR/variants.txt | tr '\n' '|' | sed 's/|$//')
INVALID_VARIANTS=$(tail -n +2 $OUT_DIR/tags.txt | awk '{print $1}' | sort -u | grep -Ev "^($INPUT_VARIANTS)$" | wc -l)
if [ $INVALID_VARIANTS -gt 0 ]; then
    echo "❌ Found variant indices not in input list"
    exit 1
fi

# Compare with R implementation
echo "➤ Running R reference implementation..."
Rscript ../../scripts/find_tags_from_vcor.R \
  $OUT_DIR/plink_tabular.vcor \
  $OUT_DIR/compressed.vars.txt \
  $OUT_DIR/variants.txt \
  $THRESHOLD \
  $OUT_DIR/tags_r.txt

# Sort both outputs for comparison (excluding headers)
sort -k1,1n -k2,2n <(tail -n +2 $OUT_DIR/tags.txt) > $OUT_DIR/tags_sorted.txt
sort -k1,1n -k2,2n <(tail -n +2 $OUT_DIR/tags_r.txt) > $OUT_DIR/tags_r_sorted.txt

# Compare
if diff -q $OUT_DIR/tags_sorted.txt $OUT_DIR/tags_r_sorted.txt > /dev/null; then
    echo "✅ C++ and R implementations match!"
else
    echo "❌ Mismatch between C++ and R implementations"
    echo "C++ output:"
    head -20 $OUT_DIR/tags_sorted.txt
    echo "R output:"
    head -20 $OUT_DIR/tags_r_sorted.txt
    echo "Diff:"
    diff $OUT_DIR/tags_sorted.txt $OUT_DIR/tags_r_sorted.txt | head -20
    exit 1
fi

# Verify at least some tags found
if [ $TAG_COUNT -gt 0 ]; then
    echo "✅ All checks passed!"
    echo "   - Tag count: $TAG_COUNT"
    echo "   - C++ vs R: identical"
fi
