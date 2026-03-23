#!/bin/bash
# Test: concat --overlapping for tabularPlink (FULL format)
#
# Strategy
# --------
# 1. Select N variants from chr20 (by position order).
# 2. Split them into two overlapping subsets:
#      chunk_a : first  A variants  (indices 0 .. A-1)
#      chunk_b : last   B variants  (indices N-B .. N-1)
#    The overlap region contains the variants shared by both subsets.
# 3. Run plink2 to generate vcor for each subset.
# 4. Compress each with plinkTabular.
# 5. Run `concat --overlapping` to merge them.
# 6. Decompress and compare the result against the full-N vcor (roundtrip accuracy).

set -e
OUT_DIR=output
mkdir -p $OUT_DIR
rm -rf $OUT_DIR/*
source ./check_plink.sh

echo -e "\n\033[1;33m ----> Test LDZip Overlapping Concat \033[0m"

BITS_LIST=(8 16 32 99)
MIN_LIST=(0.01 0.1)
N=200        # total variants
A=150        # variants in chunk_a (first A)
B=150        # variants in chunk_b (last B)
# Overlap size = A + B - N = 100

get_threshold() {
  local bits="$1"
  local min="$2"
  awk -v b="$bits" -v m="$min" '$1 == b && $2 == m {print $3}' assets/thresholds.tsv
}

# -------------------------------------------------------------------------
# Build per-chunk snplists from chr20 (skip the header line)
# -------------------------------------------------------------------------
echo "➤ Building variant subsets (N=$N, chunk_a=$A, chunk_b=$B) ..."

awk 'NR>1 {print $3}' ../../assets/g1k.chr20.pvar | head -n $N > $OUT_DIR/all.snps
head -n $A  $OUT_DIR/all.snps > $OUT_DIR/snps_a.txt
tail -n $B  $OUT_DIR/all.snps > $OUT_DIR/snps_b.txt

# -------------------------------------------------------------------------
# Generate pvar files for each chunk (needed by plinkTabular as --snp_file)
# -------------------------------------------------------------------------
for chunk in a b; do
    ${PLINK2} \
        --pfile ../../assets/g1k.chr20 \
        --extract $OUT_DIR/snps_${chunk}.txt \
        --make-pgen \
        --out $OUT_DIR/chunk_${chunk} > /dev/null
done

for min in "${MIN_LIST[@]}"; do
    echo "➤ Generating vcor files (min=$min) ..."
    min_r2=$(echo "$min * $min" | bc -l)

    for chunk in a b; do
        ${PLINK2} \
            --pfile $OUT_DIR/chunk_${chunk} \
            --ld-window-kb 1000 \
            --ld-window-r2 $min_r2 \
            --r-unphased ref-based cols=id,ref,alt \
            --out $OUT_DIR/chunk_${chunk} > /dev/null
    done

    # Also generate the full N-variant reference vcor for roundtrip comparison
    ${PLINK2} \
        --pfile ../../assets/g1k.chr20 \
        --extract $OUT_DIR/all.snps \
        --ld-window-kb 1000 \
        --ld-window-r2 $min_r2 \
        --r-unphased ref-based cols=id,ref,alt \
        --out $OUT_DIR/full_ref > /dev/null

    for bits in "${BITS_LIST[@]}"; do
        echo
        echo -e "\033[1;33m>>> Testing bits=$bits, min=$min \033[0m"

        echo "➤ Compressing chunks ..."
        for chunk in a b; do
            ../bin/ldzip compress plinkTabular \
                --ld_file   $OUT_DIR/chunk_${chunk}.vcor \
                --snp_file  $OUT_DIR/chunk_${chunk}.pvar \
                --output_prefix $OUT_DIR/compressed_${chunk} \
                --min  $min \
                --bits $bits \
                --min_col UNPHASED_R > /dev/null
        done

        echo "➤ Overlapping concat ..."
        ../bin/ldzip concat \
            --inputs $OUT_DIR/compressed_a $OUT_DIR/compressed_b \
            --output_prefix $OUT_DIR/overlap_concat \
            --overlapping > /dev/null

        # Verify variant count in the merged matrix equals N (no duplicates)
        merged_n=$(python3 -c "
import json, sys
with open('$OUT_DIR/overlap_concat.meta.json') as f:
    m = json.load(f)
print(m['dim']['rows'])
")
        echo "➤ Checking merged variant count: $merged_n (expected $N)"
        if [ "$merged_n" -ne "$N" ]; then
            echo "❌ Variant count mismatch: got $merged_n, expected $N"
            exit 1
        fi
        echo "✅ Variant count correct"

        echo "➤ Decompressing ..."
        ../bin/ldzip decompress \
            --input_prefix $OUT_DIR/overlap_concat \
            --output_prefix $OUT_DIR/overlap_full \
            --type tabular > /dev/null

        echo "➤ Comparing roundtrip against full reference ..."
        Rscript ../../scripts/corr.R \
            -orig $OUT_DIR/full_ref.vcor \
            -new  $OUT_DIR/overlap_full.vcor > /dev/null

        max_diff=$(jq .max_diff "$OUT_DIR/overlap_full_stats.json")
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

echo -e "\n\033[1;32m ✅ All overlapping concat tests passed \033[0m"
