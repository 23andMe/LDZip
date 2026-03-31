#!/bin/bash
set -e
OUT_DIR=output
rm -rf $OUT_DIR/*
N=1000
source ./check_plink.sh

echo -e "\n\033[1;33m ----> Test LDZip Tabular Concat \033[0m"

BITS_LIST=(8 16 32 99)
MIN_LIST=(0.01 0.1)

get_threshold() {
  local bits="$1"
  local min="$2"
  awk -v b="$bits" -v m="$min" '$1 == b && $2 == m {print $3}' assets/thresholds.tsv
}

# overlapping chunks on chr20 (overlaps ~45kb each, ld-window-kb 40)
# chunk0 ends at 155000, chunk2 starts at 165000 => no non-adjacent overlap
CHUNKS=(
  "60343 155000"
  "110000 210000"
  "165000 227860"
)

FULL_START=60343
FULL_END=227860

for min in "${MIN_LIST[@]}"; do
    echo "➤ Generating random plink matrices (chr20 chunks)..."
    min_r2=$(echo "$min * $min" | bc -l)

    # Full-range reference for roundtrip comparison
    ${PLINK2} \
    --pfile ../../assets/g1k.chr20 \
    --chr 20 \
    --from-bp $FULL_START \
    --to-bp $FULL_END \
    --ld-window-kb 40 \
    --ld-window-r2 $min_r2 \
    --r-unphased ref-based cols=id,ref,alt \
    --out $OUT_DIR/chr20_full > /dev/null

    awk 'NR==1 || ($2 >= '"$FULL_START"' && $2 <= '"$FULL_END"')' \
        ../../assets/g1k.chr20.pvar > $OUT_DIR/chr20_full.pvar

    idx=0
    for chunk in "${CHUNKS[@]}"; do
        start=$(echo $chunk | awk '{print $1}')
        end=$(echo $chunk | awk '{print $2}')

        ${PLINK2} \
        --pfile ../../assets/g1k.chr20 \
        --chr 20 \
        --from-bp $start \
        --to-bp $end \
        --ld-window-kb 40 \
        --ld-window-r2 $min_r2 \
        --r-unphased ref-based cols=id,ref,alt \
        --out $OUT_DIR/chr20_chunk${idx} > /dev/null

        awk 'NR==1 || ($2 >= '"$start"' && $2 <= '"$end"')' \
            ../../assets/g1k.chr20.pvar > $OUT_DIR/chr20_chunk${idx}.pvar

        idx=$((idx+1))
    done

    for bits in "${BITS_LIST[@]}"; do
        echo
        echo -e "\033[1;33m>>> Testing bits=$bits, min=$min \033[0m"

        echo "➤ Compressing..."
        idx=0
        for chunk in "${CHUNKS[@]}"; do
            ../bin/ldzip compress plinkTabular \
            --ld_file $OUT_DIR/chr20_chunk${idx}.vcor \
            --snp_file $OUT_DIR/chr20_chunk${idx}.pvar \
            --output_prefix $OUT_DIR/compressed.chunk${idx} \
            --min $min \
            --bits $bits \
            --min_col UNPHASED_R > /dev/null

            idx=$((idx+1))
        done

        echo "➤ Concating..."
        ../bin/ldzip concat \
        --inputs \
            $OUT_DIR/compressed.chunk0 \
            $OUT_DIR/compressed.chunk1 \
            $OUT_DIR/compressed.chunk2 \
        --output_prefix $OUT_DIR/concat \
        --overlapping > /dev/null

        echo "➤ Decompressing ..."
        ../bin/ldzip decompress \
        --input_prefix $OUT_DIR/concat \
        --output_prefix $OUT_DIR/concat_full \
        --type tabular > /dev/null

        echo "➤ Comparing roundtrip..."
        Rscript ../../scripts/corr.R \
            -orig $OUT_DIR/chr20_full.vcor \
            -new $OUT_DIR/concat_full.vcor > /dev/null

        max_diff=$(jq .max_diff "$OUT_DIR/concat_full_stats.json")
        max_allowed=$(get_threshold "$bits" "$min")

        echo "➤ Check max_diff[$max_diff] <= threshold[$max_allowed]"
        if awk "BEGIN {exit !($max_diff <= $max_allowed)}"; then
            echo "✅ Within threshold"
        else
            echo "❌ Exceeds threshold"
            exit 1
        fi
    done
done

echo -e "\n\033[1;32m✅ Tabular concat overlap test passed! \033[0m"