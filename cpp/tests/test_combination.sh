#!/bin/bash
set -euo pipefail
OUT_DIR=output
N=1000

echo -e "\n\033[1;33m ----> Test Lossy Compression \033[0m"

echo "➤ Generating random matrix ..."
Rscript ../../scripts/generate_random_matrix.R $N $OUT_DIR/test_in

BITS_LIST=(8 16 32 99)
MIN_LIST=(0.0 0.0001 0.001 0.01)

get_threshold() {
  local bits="$1"
  local min="$2"
  awk -v b="$bits" -v m="$min" '$1 == b && $2 == m {print $3}' assets/thresholds.tsv
}

for bits in "${BITS_LIST[@]}"; do
  for min in "${MIN_LIST[@]}"; do
    echo
    echo -e "\033[1;33m>>> Testing bits=$bits, min=$min \033[0m"

    PREFIX="${OUT_DIR}/compressed_b${bits}_min${min}"
    OUT_STATS="${OUT_DIR}/test_out_stats.json"

    echo "➤ Compressing..."
    ../bin/ldzip compress plinkSquare \
      --ld_file $OUT_DIR/test_in.bin \
      --snp_file $OUT_DIR/test_in.bin.vars \
      --output_prefix $PREFIX \
      --min $min \
      --bits $bits > /dev/null

    echo "➤ Decompressing..."
    ../bin/ldzip decompress \
      --input_prefix $PREFIX \
      --output_prefix $OUT_DIR/test_out > /dev/null

    echo "➤ Comparing roundtrip..."
    Rscript ../../scripts/corr.R \
      -orig $OUT_DIR/test_in.bin \
      -new $OUT_DIR/test_out.bin > /dev/null

    max_diff=$(jq .max_diff "$OUT_STATS")
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
