#!/bin/bash
# Test script for UKBB SuSiE analysis tutorial

set -e

# Create directory
mkdir -p test_output/ukbb_susie

# Download GWAS summary statistics
wget -P test_output/ukbb_susie https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/additive-tsvs/20002_1111.gwas.imputed_v3.both_sexes.tsv.bgz

# Download variant annotations
wget -P test_output/ukbb_susie https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/annotations/variants.tsv.bgz

# Define region (GATA3 locus)
CHR=10
START=7750000
END=8350000
REGION="chr${CHR}_${START}_${END}"

mkdir -p test_output/ukbb_susie/regions

# Extract GWAS summary statistics
bgzip -dc test_output/ukbb_susie/20002_1111.gwas.imputed_v3.both_sexes.tsv.bgz | \
    awk -F'\t' -v chr="$CHR" -v start="$START" -v end="$END" '
      NR==1 {print; next}
      {
        split($1, a, ":")
        if (a[1] == chr && a[2] >= start && a[2] <= end) print
      }' > test_output/ukbb_susie/regions/${REGION}_gwas.tsv

# Extract variant annotations
tabix -s 2 -b 3 -e 3 -S 1 test_output/ukbb_susie/variants.tsv.bgz
bgzip -dc test_output/ukbb_susie/variants.tsv.bgz | head -1 > test_output/ukbb_susie/regions/${REGION}_vars.tsv
tabix test_output/ukbb_susie/variants.tsv.bgz ${CHR}:${START}-${END} >> test_output/ukbb_susie/regions/${REGION}_vars.tsv

# Run analysis
# Note: This requires ukbb_ld/european_std to exist.
# Generate it first by running test_ukbb.sh or following the UKBB tutorial.
# See prerequisites in docs/tutorials/ukbb-susie-analysis.md
Rscript -e "
library(data.table)
library(susieR)
library(LDZipMatrix)

chr <- 10; start <- 7750000; end <- 8350000
region_name <- sprintf('chr%d_%d_%d', chr, start, end)

gwas_file <- sprintf('test_output/ukbb_susie/regions/%s_gwas.tsv', region_name)
gwas <- fread(gwas_file)
print(sprintf('Loaded %d GWAS variants', nrow(gwas)))

vars_file <- sprintf('test_output/ukbb_susie/regions/%s_vars.tsv', region_name)
vars <- fread(vars_file)
print(sprintf('Loaded %d variant annotations', nrow(vars)))

gwas <- merge(gwas, vars[, .(variant, rsid)], by = 'variant')

region_query <- sprintf('%d:%d-%d', chr, start, end)
ld_pointer <- LDZipMatrix('ukbb_ld/european_std')
ld_mat <- fetchLD(ld_pointer, region_query, region_query)
print(sprintf('Loaded %d x %d LD matrix', nrow(ld_mat), ncol(ld_mat)))

common_rsids <- intersect(gwas\$rsid, rownames(ld_mat))
print(sprintf('%d variants in common', length(common_rsids)))

gwas <- gwas[match(common_rsids, gwas\$rsid)]
ld_mat <- ld_mat[common_rsids, common_rsids]
print(sprintf('Final: %d variants', nrow(gwas)))

z <- gwas\$beta / gwas\$se
res <- susie_rss(z, ld_mat, n = median(gwas\$n_complete_samples), L = 10)

gwas[, pip := res\$pip]
cs <- susie_get_cs(res, coverage = 0.95)

print(sprintf('Found %d credible sets', length(cs\$cs)))
if (length(cs\$cs) > 0) {
  for (i in seq_along(cs\$cs)) {
    top_idx <- cs\$cs[[i]][which.max(res\$pip[cs\$cs[[i]]])]
    print(sprintf('Credible Set %d: %d variants, top variant %s (PIP=%.3f, p=%.2e)',
                  i, length(cs\$cs[[i]]), gwas\$rsid[top_idx], res\$pip[top_idx], gwas\$pval[top_idx]))
  }
}
"

echo "UKBB SuSiE analysis test completed successfully!"
