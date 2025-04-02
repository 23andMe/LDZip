#!/bin/bash
set -e
OUT_DIR=output
rm -rf $OUT_DIR/*
N=1000
source ./check_plink.sh

echo -e "\n\033[1;33m ----> Test LDZip Square Concat \033[0m"

BITS_LIST=(8 16 32 99)
MIN_LIST=(0.0 0.01 0.1)

get_threshold() {
  local bits="$1"
  local min="$2"
  awk -v b="$bits" -v m="$min" '$1 == b && $2 == m {print $3}' assets/thresholds.tsv
}

echo "➤ Generating random plink matrices ..."
for chr in 20 21 22; do
    ${PLINK2} \
    --pfile ../../assets/g1k.chr${chr} \
    --geno 0 \
    --write-snplist \
    --out output/temp_filtered > /dev/null 
    sort -R output/temp_filtered.snplist | head -n "$N" > output/top.snps
    ${PLINK2} \
    --pfile ../../assets/g1k.chr${chr} \
    --extract  output/top.snps \
    --snps-only \
    --r-unphased bin4 ref-based \
    --out $OUT_DIR/chr${chr} > /dev/null
done

for min in "${MIN_LIST[@]}"; do
    for bits in "${BITS_LIST[@]}"; do
    echo
    echo -e "\033[1;33m>>> Testing bits=$bits, min=$min \033[0m"

    TYPE=PHASED_R
    echo "➤ Compressing using type [$TYPE]..."
    for chr in 20 21 22; do
        ../bin/ldzip compress plinkSquare \
        --ld_file $OUT_DIR/chr${chr}.unphased.vcor1.bin \
        --snp_file $OUT_DIR/chr${chr}.unphased.vcor1.bin.vars \
        --output_prefix $OUT_DIR/compressed.chr${chr} \
        --min $min \
        --bits $bits \
        --type $TYPE > /dev/null
    done

    echo "➤ Concating..."
    ../bin/ldzip concat \
    --inputs $OUT_DIR/compressed.chr20 $OUT_DIR/compressed.chr21 $OUT_DIR/compressed.chr22 \
    --output_prefix $OUT_DIR/concat > /dev/null

    echo "➤ Decompressing ..."
    ../bin/ldzip decompress \
    --input_prefix $OUT_DIR/concat \
    --output_prefix $OUT_DIR/concat_full > /dev/null

    echo "➤ Comparing roundtrip..."
    mv $OUT_DIR/concat_full.$TYPE.bin $OUT_DIR/concat_full.bin
    Rscript ../../scripts/corr.R \
        -orig $OUT_DIR/chr20.unphased.vcor1.bin $OUT_DIR/chr21.unphased.vcor1.bin $OUT_DIR/chr22.unphased.vcor1.bin\
        -new $OUT_DIR/concat_full.bin > /dev/null

    max_diff=$(jq .max_diff "$OUT_DIR/concat_full_stats.json")
    max_allowed=$(get_threshold "$bits" "$min")

    echo "➤ Check max_diff[$max_diff] <= threshold[$max_allowed]"
    if awk "BEGIN {exit !($max_diff <= $max_allowed)}"; then
        echo "✅ Within threshold"
    else
        echo " ❌ Exceeds threshold"
        exit 1
    fi
  done
done
