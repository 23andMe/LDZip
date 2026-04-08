#!/bin/bash
set -e

OUT_DIR=output
rm -rf $OUT_DIR/*
mkdir -p $OUT_DIR

N=1000000
M=500
OVERLAP=150000

echo -e "\n\033[1;33m ----> Test Overflow Concat with Overlapping (2 chunks) \033[0m"

# Generate base matrix with many entries
echo "➤ Generating random matrix ..."
Rscript ../../scripts/generate_large_random_tabular.R ${N} $OUT_DIR/test_base ${M} > /dev/null

# Define 2 chunks with large overlap
# Chunk 1: rs1 to rs500000
# Chunk 2: rs350001 to rs1000000 (150K overlap: rs350001-rs500000)
CHUNK1_START=1
CHUNK1_END=500000
CHUNK2_START=$((CHUNK1_END - OVERLAP + 1))
CHUNK2_END=1000000

# Filter: keep only entries where both IDs are entirely within one chunk
# This removes entries that span across exclusive regions
awk -v c1s=$CHUNK1_START -v c1e=$CHUNK1_END -v c2s=$CHUNK2_START -v c2e=$CHUNK2_END '
NR==1 {print; next}
{
  id_a = substr($1, 3) + 0
  id_b = substr($4, 3) + 0

  # Check if both in chunk1
  in_chunk1 = (id_a >= c1s && id_a <= c1e && id_b >= c1s && id_b <= c1e)

  # Check if both in chunk2
  in_chunk2 = (id_a >= c2s && id_a <= c2e && id_b >= c2s && id_b <= c2e)

  if (in_chunk1 || in_chunk2) print
}
' $OUT_DIR/test_base.vcor > $OUT_DIR/test_full.vcor

cp $OUT_DIR/test_base.vars $OUT_DIR/test_full.vars

# Extract chunks
echo "➤ Extracting chunks..."
for chunk_num in 1 2; do
  if [ $chunk_num -eq 1 ]; then
    start=$CHUNK1_START
    end=$CHUNK1_END
  else
    start=$CHUNK2_START
    end=$CHUNK2_END
  fi

  prefix="$OUT_DIR/test_chunk${chunk_num}"

  # Extract vars
  awk -v s=$start -v e=$end 'NR==1 || ($3 ~ /^rs[0-9]+$/ && substr($3,3)+0 >= s && substr($3,3)+0 <= e)' \
    $OUT_DIR/test_full.vars > ${prefix}.vars

  # Extract vcor
  awk -v s=$start -v e=$end '
    NR==1 {print; next}
    {
      id_a = substr($1, 3) + 0
      id_b = substr($4, 3) + 0
      if (id_a >= s && id_a <= e && id_b >= s && id_b <= e) print
    }
  ' $OUT_DIR/test_full.vcor > ${prefix}.vcor
done
echo "➤ Compressing chunks..."
for k in 1 2; do
  ../bin/ldzip compress plinkTabular \
    --ld_file $OUT_DIR/test_chunk${k}.vcor \
    --snp_file $OUT_DIR/test_chunk${k}.vars \
    --output_prefix $OUT_DIR/compressed_${k} \
    --min 0.0 \
    --bits 99 \
    --min_col UNPHASED_R > /dev/null
done


echo "➤ Concating chunks..."
../bin/ldzip concat \
  --inputs $OUT_DIR/compressed_1 $OUT_DIR/compressed_2 \
  --output_prefix $OUT_DIR/concat \
  --overlapping > /dev/null

echo "➤ Compressing full matrix..."
../bin/ldzip compress plinkTabular \
  --ld_file $OUT_DIR/test_full.vcor \
  --snp_file $OUT_DIR/test_full.vars \
  --output_prefix $OUT_DIR/compressed_full \
  --min 0.0 \
  --bits 99 \
  --min_col UNPHASED_R > /dev/null


echo "➤ Decompressing concat result..."
../bin/ldzip decompress \
  --input_prefix $OUT_DIR/concat \
  --output_prefix $OUT_DIR/concat_out \
  --type tabular > /dev/null

echo "➤ Decompressing full matrix..."
../bin/ldzip decompress \
  --input_prefix $OUT_DIR/compressed_full \
  --output_prefix $OUT_DIR/full_out \
  --type tabular > /dev/null

echo "➤ Comparing concat with full..."
Rscript ../../scripts/corr.R \
    -orig $OUT_DIR/full_out.vcor \
    -new $OUT_DIR/concat_out.vcor > /dev/null
./test_identical.sh $OUT_DIR/concat_out_stats.json

echo -e "\n\033[1;32m✅ Overflow concat overlap test passed! \033[0m"
