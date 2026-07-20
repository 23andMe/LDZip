#!/bin/bash
# Test script for UKBB tutorial quick run

set -e

# Download single chunk for testing
mkdir -p test_output/ukbb/data

S3_BUCKET="s3://broad-alkesgroup-ukbb-ld/UKBB_LD"
chunk="chr20_1_3000001"

aws s3 cp --no-sign-request "${S3_BUCKET}/${chunk}.npz" test_output/ukbb/data/
aws s3 cp --no-sign-request "${S3_BUCKET}/${chunk}.gz" test_output/ukbb/data/

# Create YAML config
cat > test_output/ukbb/ukbb.yaml <<EOF
npz_template: '\${launchDir}/test_output/ukbb/data/chr{CHR}_{CHUNK}.npz'
outdir: 'test_output/ukbb/output'
prefix: 'european_ukbb'
chroms: '20'
npz_ld_type: 'UNPHASED_R'
concat_pairwise: true
min: 0.1
bits: 8
EOF

# Export ldzip path if needed
export LDZIP=$(pwd)/../../cpp/bin/ldzip

# Run pipeline
nextflow run ../../pipelines/wholeGenomeLD/main.nf -params-file test_output/ukbb/ukbb.yaml -resume

# Verify output
Rscript -e "
library(LDZipMatrix)
ld <- LDZipMatrix('test_output/ukbb/output/whole_genome/european_ukbb')
result <- fetchLD(ld, 'rs995008', 'rs4813467')
cat('LD between rs995008 and rs4813467:', result, '\n')
"

echo "UKBB test completed successfully!"
