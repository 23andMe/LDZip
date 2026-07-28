# UK Biobank SuSiE Analysis with LDZip-compressed LD Matrices

This guide demonstrates how to perform fine-mapping on UK Biobank GWAS data (from [Neale Lab](http://www.nealelab.is/uk-biobank)) using [SuSiE](https://github.com/stephenslab/susieR) with LDZip-compressed LD matrices (from [Alkes Price Lab](https://labs.icahn.mssm.edu/minervalab/resources/data-ark/ukbb_ld/)).

**Prerequisites**:
- R packages (`data.table`, `susieR`, `LDZipMatrix`)
- Whole-genome LD file for UKBB at `ukbb_ld/european_std` ([see UKBB LD tutorial](ukbb-tutorial.md))
- Command-line tools (`bgzip`, `tabix`)

## Workflow Overview

This example fine-maps **Asthma (phenotype 20002_1111)** at the ***GATA3* locus (chr10:7.75-8.35 Mb, GRCh38)**.

1. Download GWAS summary statistics and variant annotations
2. Extract region from whole-genome files  
3. Run SuSiE analysis in R

## Step 1: Download UK Biobank Data for given Phenotype

```bash
# Create directory
mkdir -p ukbb

# Download GWAS summary statistics (example: Asthma, phenotype 20002_1111)
wget -P ukbb https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/additive-tsvs/20002_1111.gwas.imputed_v3.both_sexes.tsv.bgz

# Download variant annotations with rsIDs
wget -P ukbb https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/annotations/variants.tsv.bgz
```

## Step 2: Extract Region

```bash
# Define region (GATA3 locus: chr10:7.75-8.35 Mb, GRCh38)
CHR=10
START=7750000
END=8350000
REGION="chr${CHR}_${START}_${END}"

# Create output directory
mkdir -p ukbb/regions

# Extract GWAS summary statistics using awk
bgzip -dc ukbb/20002_1111.gwas.imputed_v3.both_sexes.tsv.bgz | \
    awk -F'\t' -v chr="$CHR" -v start="$START" -v end="$END" '
      NR==1 {print; next}
      {
        split($1, a, ":")
        if (a[1] == chr && a[2] >= start && a[2] <= end) print
      }' > ukbb/regions/${REGION}_gwas.tsv

# Extract variant annotations using tabix (will create index if needed)
tabix -s 2 -b 3 -e 3 -S 1 ukbb/variants.tsv.bgz
bgzip -dc ukbb/variants.tsv.bgz | head -1 > ukbb/regions/${REGION}_vars.tsv
tabix ukbb/variants.tsv.bgz ${CHR}:${START}-${END} >> ukbb/regions/${REGION}_vars.tsv
```

## Step 3: Run Analysis in R

The R script loads the pre-extracted files from disk (fast) and performs the analysis:

```r
#!/usr/bin/env Rscript
# UK Biobank SuSiE Fine-mapping with LDZip
# Note: Run region extraction first to subset data

library(data.table)
library(susieR)
library(LDZipMatrix)

# Define region (GATA3 locus: chr10:7.75-8.35 Mb, GRCh38)
chr <- 10; start <- 7750000; end <- 8350000
region_name <- sprintf("chr%d_%d_%d", chr, start, end)

# Load pre-extracted GWAS data
gwas_file <- sprintf("ukbb/regions/%s_gwas.tsv", region_name)
gwas <- fread(gwas_file)
print(sprintf("Loaded %d GWAS variants", nrow(gwas)))

# Load pre-extracted variant annotations
vars_file <- sprintf("ukbb/regions/%s_vars.tsv", region_name)
vars <- fread(vars_file)
print(sprintf("Loaded %d variant annotations", nrow(vars)))

# Merge to decorate with rsIDs
gwas <- merge(gwas, vars[, .(variant, rsid, ref, alt)], by = "variant")

# Load LD matrix from LDZip
region_query <- sprintf("%d:%d-%d", chr, start, end)
ld_pointer <- LDZipMatrix("ukbb_ld/european_std")
ld_mat <- fetchLD(ld_pointer, region_query, region_query)
ld_vars = fetchVariants(ld_pointer, region_query) 
print(sprintf("Loaded %d x %d LD matrix", nrow(ld_mat), ncol(ld_mat)))

# Intersect and subset data
common_rsids <- intersect(gwas$rsid, rownames(ld_mat))
print(sprintf("%d variants in common", length(common_rsids)))

gwas <- gwas[match(common_rsids, gwas$rsid)]
ld_mat <- ld_mat[common_rsids, common_rsids]
ld_vars = ld_vars[match(common_rsids, ld_vars$ID),]
print(sprintf("Final: %d variants", nrow(gwas)))

# Harmonize alleles: flip effect size if GWAS alt != LD matrix ALT
# Note: In this case, GWAS and LD matrix use the same allele coding, so no flipping occurs.
# However, this check is good practice to always include when working with different data sources.
flip <- gwas$alt != ld_vars$ALT
print(sprintf("Flipping %d/%d variants where GWAS alt != LD ALT", sum(flip), length(flip)))

# Run SuSiE
z <- gwas$beta / gwas$se
z[flip] <- -z[flip]
res <- susie_rss(z, ld_mat, n = median(gwas$n_complete_samples), L = 10)

# Extract credible sets
gwas[, pip := res$pip]
cs <- susie_get_cs(res, coverage = 0.95)

# Summary of results
print(sprintf("Found %d credible sets", length(cs$cs)))
if (length(cs$cs) > 0) {
  for (i in seq_along(cs$cs)) {
    top_idx <- cs$cs[[i]][which.max(res$pip[cs$cs[[i]]])]
    print(sprintf("Credible Set %d: %d variants, top variant %s (PIP=%.3f, p=%.2e)",
                  i, length(cs$cs[[i]]), gwas$rsid[top_idx], res$pip[top_idx], gwas$pval[top_idx]))
  }
}

# Expected output:
# [1] "Loaded 3704 GWAS variants"
# [1] "Loaded 3704 variant annotations"
# [1] "Loaded 4957 x 4957 LD matrix"
# [1] "3633 variants in common"
# [1] "Final: 3633 variants"
# [1] "Found 3 credible sets"
# [1] "Credible Set 1: 1 variants, top variant rs11567923 (PIP=0.961, p=1.39e-10)"
# [1] "Credible Set 2: 51 variants, top variant rs117158080 (PIP=0.080, p=1.21e-07)"
# [1] "Credible Set 3: 9 variants, top variant rs263424 (PIP=0.406, p=2.49e-05)"
```

## Troubleshooting

**If you get 'XtX is not symmetric' warning:**

This can occur due to small asymmetries in the LD matrix from quantization of overlapping regions. See the [UKBB LD tutorial](ukbb-tutorial.md) for details on matrix asymmetry. This warning should not affect the final results.

## References

- [Neale Lab UK Biobank GWAS](http://www.nealelab.is/uk-biobank)
- [Alkes Price Lab UK Biobank LD](https://labs.icahn.mssm.edu/minervalab/resources/data-ark/ukbb_ld/)
- [SuSiE (Sum of Single Effects)](https://github.com/stephenslab/susieR)
