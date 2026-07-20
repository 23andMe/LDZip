#!/bin/bash
# Test script for 1000 Genomes tutorial quick run

set -e

# Download VCF for chr20
mkdir -p test_output/g1k/data

BASE_URL="ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased"
FILE_PREFIX="CCDG_14151_B01_GRM_WGS_2020-08-05"
chr=20

wget -O test_output/g1k/data/1000g.chr${chr}.vcf.gz ${BASE_URL}/${FILE_PREFIX}_chr${chr}.filtered.shapeit2-duohmm-phased.vcf.gz
wget -O test_output/g1k/data/1000g.chr${chr}.vcf.gz.tbi ${BASE_URL}/${FILE_PREFIX}_chr${chr}.filtered.shapeit2-duohmm-phased.vcf.gz.tbi

# Download sample panel
wget -P test_output/g1k/data https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel

# Retrieve list of EUR samples
awk '$3=="EUR" {print $1}' test_output/g1k/data/integrated_call_samples_v3.20130502.ALL.panel > test_output/g1k/data/EUR.txt

# Create YAML config
cat > test_output/g1k/1000g.yaml <<EOF
vcf_template: '\${launchDir}/test_output/g1k/data/1000g.chr{CHR}'
keep: '\${launchDir}/test_output/g1k/data/EUR.txt'
outdir: 'test_output/g1k/output'
ld_command: '--r-unphased ref-based cols=id,ref,alt'
prefix: 'EUR'
chroms: '20'
ld_window_kb: 1000
ld_window_r2: 0.01
min_col: 'UNPHASED_R'
ld_threads: 1
chunk_size_kb: 20000
overlap_size_kb: 1000
EOF

# Export ldzip path if needed
export LDZIP=$(pwd)/../../cpp/bin/ldzip

# Run pipeline
nextflow run ../../pipelines/wholeGenomeLD/main.nf -params-file test_output/g1k/1000g.yaml -resume

# Verify output
Rscript -e "
library(LDZipMatrix)
ld <- LDZipMatrix('test_output/g1k/output/whole_genome/EUR')
result <- fetchLD(ld, '20:64331475:C:T', '20:64333832:C:A')
cat('LD between 20:64331475:C:T and 20:64333832:C:A:', result, '\n')
"

echo "1000 Genomes test completed successfully!"
