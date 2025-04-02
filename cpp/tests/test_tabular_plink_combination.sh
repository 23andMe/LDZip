#!/bin/bash
set -e
OUT_DIR=output
rm -rf $OUT_DIR/*
source ./check_plink.sh

echo -e "\n\033[1;33m ----> Test PLINK Tabular Compression \033[0m"

BITS_LIST=(8 16 32 99)
MIN_LIST=(0.01 0.1)

get_threshold() {
  local bits="$1"
  local min="$2"
  awk -v b="$bits" -v m="$min" '$1 == b && $2 == m {print $3}' assets/thresholds.tsv
}

for min in "${MIN_LIST[@]}"; do
  echo "➤ Generating random plink matrix ..."
  min_r2=$(echo "$min * $min" | bc -l)
  ${PLINK2} \
    --pfile ../../assets/g1k \
    --ld-window-kb 100 \
    --ld-window-r2 $min_r2 \
    --r-unphased ref-based cols=id,ref,alt \
    --out $OUT_DIR/plink_tabular > /dev/null

 
  for bits in "${BITS_LIST[@]}"; do
    echo
    echo -e "\033[1;33m>>> Testing bits=$bits, min=$min \033[0m"

    PREFIX="${OUT_DIR}/compressed_b${bits}_min${min}"
    OUT_STATS="${OUT_DIR}/test_out_stats.json"

    echo "➤ Compressing..."
    ../bin/ldzip compress plinkTabular \
      --ld_file $OUT_DIR/plink_tabular.vcor \
      --snp_file ../../assets/g1k.pvar \
      --output_prefix $PREFIX \
      --bits $bits \
      --min $min \
      --min_col UNPHASED_R > /dev/null


    echo "➤ Decompressing..."
    ../bin/ldzip decompress \
      --input_prefix $PREFIX \
      --output_prefix $OUT_DIR/test_out \
      --type tabular  > /dev/null

    echo "➤ Comparing roundtrip..."
    Rscript ../../scripts/corr.R \
        -orig $OUT_DIR/plink_tabular.vcor \
        -new $OUT_DIR/test_out.vcor > /dev/null

    max_diff=$(jq .max_diff "$OUT_STATS")
    max_allowed=$(get_threshold "$bits" "0")

    echo "➤ Check max_diff[$max_diff] <= threshold[$max_allowed]"
    if awk "BEGIN {exit !($max_diff <= $max_allowed)}"; then
        echo "✅ Within threshold"
    else
        echo " ❌ Exceeds threshold"
        exit 1
    fi


  done
done

