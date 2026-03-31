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
echo "➤ Generating base matrix (N=$N, M=$M)..."
Rscript ../../scripts/generate_large_random_tabular.R ${N} $OUT_DIR/test_base ${M} > /dev/null

# Define 2 chunks with large overlap
# Chunk 1: rs1 to rs500000
# Chunk 2: rs350001 to rs1000000 (150K overlap: rs350001-rs500000)
CHUNK1_START=1
CHUNK1_END=500000
CHUNK2_START=$((CHUNK1_END - OVERLAP + 1))
CHUNK2_END=1000000

echo "➤ Chunk boundaries:"
echo "  Chunk 1: rs${CHUNK1_START} to rs${CHUNK1_END}"
echo "  Chunk 2: rs${CHUNK2_START} to rs${CHUNK2_END}"
echo "  Overlap: rs${CHUNK2_START} to rs${CHUNK1_END} (${OVERLAP} variants)"

# Filter: keep only entries where both IDs are entirely within one chunk
# This removes entries that span across exclusive regions
echo "➤ Filtering to remove spanning entries..."
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

ORIG_COUNT=$(tail -n +2 $OUT_DIR/test_base.vcor | wc -l | tr -d ' ')
FILT_COUNT=$(tail -n +2 $OUT_DIR/test_full.vcor | wc -l | tr -d ' ')
echo "  Original entries: $ORIG_COUNT"
echo "  Filtered entries: $FILT_COUNT (removed $(($ORIG_COUNT - $FILT_COUNT)) spanning entries)"

# Analyze distribution across regions
echo "➤ Analyzing overflow distribution..."
awk -v c1s=$CHUNK1_START -v c1e=$CHUNK1_END -v c2s=$CHUNK2_START -v c2e=$CHUNK2_END '
NR==1 {next}
{
  id_a = substr($1, 3) + 0
  id_b = substr($4, 3) + 0
  delta = id_b - id_a

  # Exclusive to chunk1: both < overlap_start
  if (id_a < c2s && id_b < c2s) {
    chunk1_excl++
    if (delta > 32767) c1_overflow++
  }
  # Entirely in overlap
  else if (id_a >= c2s && id_a <= c1e && id_b >= c2s && id_b <= c1e) {
    overlap++
    if (delta > 32767) ov_overflow++
  }
  # Exclusive to chunk2: both > chunk1_end
  else if (id_a > c1e && id_b > c1e) {
    chunk2_excl++
    if (delta > 32767) c2_overflow++
  }
}
END {
  print "  Chunk1 exclusive:", chunk1_excl, "entries (" c1_overflow " overflow)"
  print "  Overlap region:", overlap, "entries (" ov_overflow " overflow)"
  print "  Chunk2 exclusive:", chunk2_excl, "entries (" c2_overflow " overflow)"
}
' $OUT_DIR/test_full.vcor

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

  vcor_count=$(tail -n +2 ${prefix}.vcor | wc -l | tr -d ' ')
  vars_count=$(tail -n +2 ${prefix}.vars | wc -l | tr -d ' ')
  echo "  Chunk ${chunk_num}: ${vars_count} variants, ${vcor_count} LD entries"
done

# Compress chunks
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

# Compress full matrix for comparison
echo "➤ Compressing full matrix..."
../bin/ldzip compress plinkTabular \
  --ld_file $OUT_DIR/test_full.vcor \
  --snp_file $OUT_DIR/test_full.vars \
  --output_prefix $OUT_DIR/compressed_full \
  --min 0.0 \
  --bits 99 \
  --min_col UNPHASED_R > /dev/null

echo "➤ Concatenating with overlapping..."
../bin/ldzip concat \
  --inputs $OUT_DIR/compressed_1 $OUT_DIR/compressed_2 \
  --output_prefix $OUT_DIR/concat \
  --overlapping > /dev/null

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

# Verify overflow handling - check that COO files exist and are non-empty
echo "➤ Verifying overflow (COO) files..."
ALL_COO_VALID=true

for k in 1 2; do
  if [ ! -f $OUT_DIR/compressed_${k}.io.bin ]; then
    echo "  ✗ ERROR: Chunk ${k} missing .io.bin file"
    ALL_COO_VALID=false
  else
    SIZE=$(stat -f%z $OUT_DIR/compressed_${k}.io.bin 2>/dev/null || stat -c%s $OUT_DIR/compressed_${k}.io.bin 2>/dev/null)
    if [ "$SIZE" -eq 0 ]; then
      echo "  ✗ ERROR: Chunk ${k} has empty .io.bin file"
      ALL_COO_VALID=false
    else
      ENTRIES=$((SIZE / 12))
      echo "  ✓ Chunk ${k}: ${ENTRIES} COO overflow entries (${SIZE} bytes)"
    fi
  fi
done

if [ ! -f $OUT_DIR/concat.io.bin ]; then
  echo "  ✗ ERROR: Concat result missing .io.bin file"
  ALL_COO_VALID=false
else
  SIZE=$(stat -f%z $OUT_DIR/concat.io.bin 2>/dev/null || stat -c%s $OUT_DIR/concat.io.bin 2>/dev/null)
  if [ "$SIZE" -eq 0 ]; then
    echo "  ✗ ERROR: Concat result has empty .io.bin file"
    ALL_COO_VALID=false
  else
    ENTRIES=$((SIZE / 12))
    echo "  ✓ Concat result: ${ENTRIES} COO overflow entries (${SIZE} bytes)"
  fi
fi

if [ "$ALL_COO_VALID" = false ]; then
  echo -e "\n\033[1;31m✗ COO overflow file validation failed! \033[0m"
  exit 1
fi

echo -e "\n\033[1;32m✅ Overflow concat overlap test passed! \033[0m"
