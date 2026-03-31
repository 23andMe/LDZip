#!/bin/bash
set -e
OUT_DIR=output
rm -rf $OUT_DIR/*
source ./check_plink.sh

min=0.01
min_r2=$(echo "$min * $min" | bc -l)

generate_and_compress() {
    local test_name=$1
    shift
    local chunks=("$@")

    idx=0
    for chunk in "${chunks[@]}"; do
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
        --out $OUT_DIR/${test_name}_chunk${idx} > /dev/null

        awk 'NR==1 || ($2 >= '"$start"' && $2 <= '"$end"')' \
            ../../assets/g1k.chr20.pvar > $OUT_DIR/${test_name}_chunk${idx}.pvar

        ../bin/ldzip compress plinkTabular \
            --ld_file $OUT_DIR/${test_name}_chunk${idx}.vcor \
            --snp_file $OUT_DIR/${test_name}_chunk${idx}.pvar \
            --output_prefix $OUT_DIR/${test_name}_c${idx} \
            --min $min --bits 16 --min_col UNPHASED_R > /dev/null

        idx=$((idx+1))
    done
}

echo -e "\n\033[1;33m ----> Test concat error handling \033[0m"

# ============================================================================
# Test 1: overlapping=T but adjacent chunks DON'T overlap → should fail
# ============================================================================
echo -e "\n\033[1;34m[Test 1] overlapping=T but chunk0 doesn't overlap with chunk1\033[0m"

CHUNKS=("60343 150000" "165000 227860")
generate_and_compress "test1" "${CHUNKS[@]}"

if ../bin/ldzip concat \
    --inputs $OUT_DIR/test1_c0 $OUT_DIR/test1_c1 \
    --output_prefix $OUT_DIR/test1_out \
    --overlapping 2>&1 | grep -q "Chunk overlap error"; then
    echo "✅ Test 1 passed"
else
    echo "❌ Test 1 failed: expected chunk overlap error"
    exit 1
fi

# ============================================================================
# Test 2: overlapping=T but chunk0 overlaps with chunk2 (non-adjacent) → should fail
# ============================================================================
echo -e "\n\033[1;34m[Test 2] overlapping=T but chunk0 overlaps with chunk2 (non-adjacent)\033[0m"

CHUNKS=("60343 170000" "110000 210000" "165000 227860")
generate_and_compress "test2" "${CHUNKS[@]}"

if ../bin/ldzip concat \
    --inputs $OUT_DIR/test2_c0 $OUT_DIR/test2_c1 $OUT_DIR/test2_c2 \
    --output_prefix $OUT_DIR/test2_out \
    --overlapping 2>&1 | grep -q "Chunk overlap error"; then
    echo "✅ Test 2 passed"
else
    echo "❌ Test 2 failed: expected chunk overlap error"
    exit 1
fi

echo -e "\n\033[1;32m✅ All concat fail tests passed! \033[0m"
