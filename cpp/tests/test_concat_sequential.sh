#!/bin/bash
set -e
OUT_DIR=output
rm -rf $OUT_DIR/*
N=1000
source ./check_plink.sh

echo -e "\n\033[1;33m ----> Test Concat Sequential (Lossless)  \033[0m"

BITS_LIST=(8 99)
MIN_LIST=(0.01 0.1)

get_threshold() {
  local bits="$1"
  local min="$2"
  awk -v b="$bits" -v m="$min" '$1 == b && $2 == m {print $3}' assets/thresholds.tsv
}

for min in "${MIN_LIST[@]}"; do
    echo "➤ Generating random plink matrices ..."
    min_r2=$(echo "$min * $min" | bc -l)
    for chr in 20 21 22; do
        ${PLINK2} \
        --pfile ../../assets/g1k.chr${chr} \
        --ld-window-kb 100 \
        --ld-window-r2 $min_r2 \
        --r-unphased ref-based cols=id,ref,alt \
        --out $OUT_DIR/chr${chr} > /dev/null
    done

    for bits in "${BITS_LIST[@]}"; do
    echo
    echo -e "\033[1;33m>>> Testing bits=$bits, min=$min \033[0m"

    echo "➤ Compressing..."
    for chr in 20 21 22; do
        ../bin/ldzip compress plinkTabular \
        --ld_file $OUT_DIR/chr${chr}.vcor \
        --snp_file ../../assets/g1k.chr${chr}.pvar \
        --output_prefix $OUT_DIR/compressed.chr${chr} \
        --min $min \
        --bits $bits \
        --min_col UNPHASED_R > /dev/null
    done

    echo "➤ Concating all together ..."
    ../bin/ldzip concat \
    --inputs $OUT_DIR/compressed.chr20 $OUT_DIR/compressed.chr21 $OUT_DIR/compressed.chr22 \
    --output_prefix $OUT_DIR/concat_all > /dev/null

    echo "➤ Concating sequentially ..."
    ../bin/ldzip concat \
    --inputs $OUT_DIR/compressed.chr20 $OUT_DIR/compressed.chr21 \
    --output_prefix $OUT_DIR/concat_1 > /dev/null
    ../bin/ldzip concat \
    --inputs $OUT_DIR/concat_1 $OUT_DIR/compressed.chr22 \
    --output_prefix $OUT_DIR/concat_sequential > /dev/null

    echo "➤ Decompressing concat_all ..."
    ../bin/ldzip decompress \
    --input_prefix $OUT_DIR/concat_all \
    --output_prefix $OUT_DIR/concat_all_full \
    --type tabular > /dev/null

    echo "➤ Decompressing concat_sequential ..."
    ../bin/ldzip decompress \
    --input_prefix $OUT_DIR/concat_sequential \
    --output_prefix $OUT_DIR/concat_sequential_full \
    --type tabular > /dev/null

    echo "➤ Comparing roundtrip..."
    Rscript ../../scripts/corr.R \
        -orig $OUT_DIR/concat_sequential_full.vcor \
        -new $OUT_DIR/concat_all_full.vcor > /dev/null

    ./test_identical.sh $OUT_DIR/concat_all_full_stats.json
  done
done
