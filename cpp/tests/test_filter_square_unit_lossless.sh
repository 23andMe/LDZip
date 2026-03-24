#!/bin/bash
set -e
OUT_DIR=output
rm -rf $OUT_DIR/*
N=1000
TYPE=PHASED_R
source ./check_plink.sh
echo -e "\n\033[1;33m ----> Test Filtering on Lossless Compression \033[0m"

echo "➤ Generating random matrix ..."
Rscript ../../scripts/generate_random_matrix.R $N $OUT_DIR/test_in
seq 0 $((N - 1)) | sort -R | head -n 500 | sort -n > "$OUT_DIR/random_keep"

echo "➤ Compressing with ldzip using type [$TYPE]..."
../bin/ldzip compress plinkSquare\
  --ld_file $OUT_DIR/test_in.bin \
  --snp_file $OUT_DIR/test_in.bin.vars \
  --output_prefix $OUT_DIR/compressed \
  --min 0.0 \
  --bits 99 \
  --type $TYPE > /dev/null

echo "➤ Filtering compressed file with random 500 indices..."
../bin/ldzip filter \
  --input_prefix $OUT_DIR/compressed \
  --keep $OUT_DIR/random_keep \
  --output_prefix $OUT_DIR/compressed_filtered  > /dev/null

echo "➤ Decompressing with ldzip..."
../bin/ldzip decompress \
  --input_prefix $OUT_DIR/compressed_filtered \
  --output_prefix $OUT_DIR/test_out  > /dev/null

echo "➤ Filtering uncompressed file..."
Rscript ../../scripts/filter.R \
    -orig $OUT_DIR/test_in.bin \
    -keep $OUT_DIR/random_keep \
    -out $OUT_DIR/test_in_filtered.bin  > /dev/null

echo "➤ Diffing original and roundtrip binary..."
diff -q <(xxd $OUT_DIR/test_in_filtered.bin) <(xxd $OUT_DIR/test_out.$TYPE.bin) \
  || (echo " ❌ Round-trip failed." && exit 1)

echo "➤ Comparing roundtrip..."
mv $OUT_DIR/test_out.$TYPE.bin $OUT_DIR/test_out.bin
Rscript ../../scripts/corr.R \
    -orig $OUT_DIR/test_in_filtered.bin \
    -new $OUT_DIR/test_out.bin  > /dev/null

./test_identical.sh $OUT_DIR/test_out_stats.json

