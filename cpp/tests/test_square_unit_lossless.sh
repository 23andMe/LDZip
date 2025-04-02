#!/bin/bash
set -e
OUT_DIR=output
rm -rf $OUT_DIR/*
N=1000
TYPE=PHASED_R
source ./check_plink.sh
echo -e "\n\033[1;33m ----> Test Lossless Compression \033[0m"

echo "➤ Generating random matrix ..."
Rscript ../../scripts/generate_random_matrix.R $N $OUT_DIR/test_in

echo "➤ Compressing with ldzip using type [$TYPE]..."
../bin/ldzip compress plinkSquare\
  --ld_file $OUT_DIR/test_in.bin \
  --snp_file $OUT_DIR/test_in.bin.vars \
  --output_prefix $OUT_DIR/compressed \
  --min 0.0 \
  --bits 99 \
  --type $TYPE > /dev/null


echo "➤ Checking compressed output files..."
for f in .x.$TYPE.bin .i.bin .io.bin .p.bin .meta.json .vars.txt; do
    if [ "$f" != ".io.bin" ] && [ ! -s "$OUT_DIR/compressed$f" ]; then
        echo "❌ Missing or empty file: $OUT_DIR/compressed$f"
        exit 1
    fi
done

echo "➤ Decompressing with ldzip..."
../bin/ldzip decompress \
  --input_prefix $OUT_DIR/compressed \
  --output_prefix $OUT_DIR/test_out


echo "➤ Checking decompressed output files..."
for f in .$TYPE.bin .bin.vars; do
    if [ ! -s "$OUT_DIR/test_out$f" ]; then
        echo "❌ Missing or empty file: $OUT_DIR/test_out$f"
        exit 1
    fi
done

echo "➤ Diffing original and roundtrip binary..."
diff -q <(xxd $OUT_DIR/test_in.bin) <(xxd $OUT_DIR/test_out.$TYPE.bin) \
  || (echo " ❌ Round-trip failed." && exit 1)

echo "➤ Comparing roundtrip..."
mv $OUT_DIR/test_out.$TYPE.bin $OUT_DIR/test_out.bin
Rscript ../../scripts/corr.R \
    -orig $OUT_DIR/test_in.bin \
    -new $OUT_DIR/test_out.bin  > /dev/null

./test_identical.sh $OUT_DIR/test_out_stats.json

