# Generating Whole-Genome LDZip Matrix from UK Biobank LD Matrices

## Overview

UKBB-LD provides summary linkage disequilibrium (LD) matrices computed from UK Biobank based on N=337K British-ancestry individuals. The LD matrices were computed by the [Alkes Price group at Harvard](https://labs.icahn.mssm.edu/minervalab/resources/data-ark/ukbb_ld/) and are publicly available in AWS S3. The LD information is stored as 2,763 3Mb-long regions spanning the entire genome in NPZ format.

If you want to use the LDZip-compressed LD matrix generated from this tutorial for SuSiE fine-mapping, see the [UKBB SuSiE analysis tutorial](ukbb-susie-analysis.md). See also the [1000 Genomes tutorial](g1k-tutorial.md) for generating LD matrices from raw VCF files.

The workflow consists of two main steps:

1. Download pre-computed NPZ LD matrices from AWS S3
2. Run a Nextflow pipeline to compress and concatenate into a whole genome LDZip matrix

A quick run option is provided to test a single chromosome chunk and ensure the workflow works correctly before scaling up to the full genome.

**Note:** The pipeline automatically handles several data-specific processing steps that are specific to the Alkes Price lab UKBB-LD data format. If you use this on other data, use it with care and verify the output:
- **Boundary variant removal**: The first variant in each 3Mb region (at position ending in 00001) is automatically removed to prevent duplication when concatenating overlapping regions.
- **Multi-allelic variant handling**: RSIDs that appear multiple times (different alleles at same position) are made unique by appending `:REF:ALT` to the variant ID.
- **Exact duplicate removal**: Variants with identical `CHR:POS:REF:ALT` are deduplicated to ensure each variant appears only once in the final matrix.
- **Diagonal correction**: Matrix diagonal values are set to 1.0 (some NPZ files have 0.5 on the diagonal, which is corrected during processing).

## Prerequisites

- [`nextflow`](https://www.nextflow.io/docs/latest/install.html)
- AWS CLI (for downloading data)
- HPC cluster (recommended for whole-genome processing)

---

## Setup: Clone and Build

Before running the pipeline, clone the repository and build the required binaries:

```bash
# Clone LDZip repository
git clone git@github.com:23andMe/LDZip.git
cd LDZip

# Build C++ binary
cd cpp
make
cd ..

# Install R package
cd R
make install
cd ..

# Return to working directory
cd ..
```

After setup, the `ldzip` binary will be at `LDZip/cpp/bin/ldzip` and the R package `LDZipMatrix` will be installed.

---

## Step 1: Download UKBB-LD NPZ Files

### Quick run

Download single chunk for testing, then proceed to Step 2.

```bash
mkdir -p data/ukbb

# Download one chr20 file (chr20: 0-3Mb)
S3_BUCKET="s3://broad-alkesgroup-ukbb-ld/UKBB_LD"
chunk="chr20_1_3000001"

aws s3 cp --no-sign-request "${S3_BUCKET}/${chunk}.npz" data/ukbb/
aws s3 cp --no-sign-request "${S3_BUCKET}/${chunk}.gz" data/ukbb/
```

### Full run

Download all chunks for all chromosomes.

```bash
mkdir -p data/ukbb

# Download all 2,763 NPZ files from AWS S3 (no AWS credentials required)
aws s3 sync --no-sign-request \
  s3://broad-alkesgroup-ukbb-ld/UKBB_LD/ \
  data/ukbb/
```

---

## Step 2: Run Nextflow LD Pipeline

### Quick run

**File: ukbb.yaml**
```yaml
npz_template: '${launchDir}/data/ukbb/chr{CHR}_{CHUNK}.npz'
outdir: 'output'
prefix: 'european_ukbb'
chroms: '20'
npz_ld_type: 'UNPHASED_R'
concat_pairwise: true
min: 0.1
bits: 8
```

If `ldzip` is NOT available in your `$PATH`, export its path:

```bash
export LDZIP=$(pwd)/LDZip/cpp/bin/ldzip
```

Run as follows:

```bash
nextflow run LDZip/pipelines/wholeGenomeLD/main.nf -params-file ukbb.yaml -resume
```

This runs locally and requires at least 16GB RAM.

Go to Step 3 to test whether it worked correctly.

### Full run

To run on all chromosomes, update the YAML:

```yaml
chroms: '1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22'
```

For whole-genome processing, you might need an HPC cluster. For example, if using SLURM:

**File: slurm.config**
```groovy
params.partition        = "example_partition"

process.executor        = "slurm"
process.cpus            = 1
process.errorStrategy   = 'retry'
process.maxRetries      = 5
process.queue           = params.partition

executor.perCpuMemAllocation = true
```

Run as follows:

```bash
nextflow run LDZip/pipelines/wholeGenomeLD/main.nf -params-file ukbb.yaml -C slurm.config -resume
```

For other HPC environments, refer to the Nextflow executor [guidelines](https://www.nextflow.io/docs/latest/executor.html).

---

## Step 3: Verify Output

After successful completion, verify the LD matrix by querying two variants:

```r
library(LDZipMatrix)

# Load the LD matrix
ld <- LDZipMatrix("output/whole_genome/european_ukbb")

# Query LD between two variants
fetchLD(ld, "rs995008", "rs4813467")
# [1] 0.5234156

# Benchmark query time
system.time(fetchLD(ld, "rs995008", "rs4813467"))
#    user  system elapsed
#   0.003   0.000   0.010
```

---

## Parameters

| Parameter        | Description                                                                 | Example / Notes |
|------------------|-----------------------------------------------------------------------------|-----------------|
| `npz_template`   | Template path to per-chromosome-chunk NPZ files. `{CHR}` and `{CHUNK}` are replaced at runtime | `${launchDir}/data/ukbb/chr{CHR}_{CHUNK}.npz` |
| `outdir`         | Directory where LD outputs will be written                                | `output` |
| `prefix`         | Prefix used for naming output files                                       | `european_ukbb` |
| `chroms`         | Comma-separated list of chromosomes to process                            | `1–22` |
| `npz_ld_type`    | LD statistic type to compress from NPZ files                              | `UNPHASED_R`, `PHASED_R`, `UNPHASED_R2`, `DPRIME` |
| `concat_pairwise`| Use pairwise concatenation for overlapping regions (required for UKBB-LD) | `true` |
| `min`            | Minimum absolute LD value to store (values below are set to zero)        | `0.1` |
| `bits`           | Number of bits for quantization (8, 16, 32, or 99 for raw float)         | `8` |
| `stage_chunk`    | *(Optional)* Stage intermediate chunk files to outdir (default: `false`)  | `true` |
| `stage_chr`      | *(Optional)* Stage per-chromosome files to outdir (default: `false`)      | `true` |
