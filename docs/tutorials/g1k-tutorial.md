# Generating Whole-Genome LDZip Matrix from 1000 Genomes 

## Overview

This tutorial demonstrates generating whole-genome linkage disequilibrium (LD) matrices from 1000 Genomes Project (1000G) phased VCF data. The same workflow can be applied to other large-scale datasets such as UK Biobank or TOPMed.

The workflow consists of two main steps:

1. Download phased VCFs for each chromosome  
2. Run a Nextflow pipeline to compute a whole genome LDZip matrix

A quick run option is provided to test a subset of chromosomes and ensure the workflow works correctly before scaling up to the full genome.

**See also:** [UK Biobank tutorial](ukbb-tutorial.md) for generating LD matrices from pre-computed NPZ files.

## Prerequisites

- [`nextflow`](https://www.nextflow.io/docs/latest/install.html)
- [`plink2`](https://www.cog-genomics.org/plink/2.0/)
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

## Step 1: Download 1000G Phased VCFs

### Quick run

Download VCF for chromosome 20 for testing, then proceed to Step 2.

```bash
mkdir -p data

# Download VCF for chr20
BASE_URL="ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased"
FILE_PREFIX="CCDG_14151_B01_GRM_WGS_2020-08-05"
chr=20

wget -O data/1000g.chr${chr}.vcf.gz ${BASE_URL}/${FILE_PREFIX}_chr${chr}.filtered.shapeit2-duohmm-phased.vcf.gz
wget -O data/1000g.chr${chr}.vcf.gz.tbi ${BASE_URL}/${FILE_PREFIX}_chr${chr}.filtered.shapeit2-duohmm-phased.vcf.gz.tbi

# Download sample panel (2504 samples, phase3)
wget -P data https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel
```

### Full run

Download VCFs for all chromosomes.

```bash
mkdir -p data

# Download VCFs for all chromosomes
BASE_URL="ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased"
FILE_PREFIX="CCDG_14151_B01_GRM_WGS_2020-08-05"

for chr in {1..22}; do
  wget -O data/1000g.chr${chr}.vcf.gz ${BASE_URL}/${FILE_PREFIX}_chr${chr}.filtered.shapeit2-duohmm-phased.vcf.gz
  wget -O data/1000g.chr${chr}.vcf.gz.tbi ${BASE_URL}/${FILE_PREFIX}_chr${chr}.filtered.shapeit2-duohmm-phased.vcf.gz.tbi
done

# Download sample panel (2504 samples, phase3)
wget -P data https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel
```


---

## Step 2: Run Nextflow LD Pipeline for European Samples

Retrieve list of EUR samples:

```bash
awk '$3=="EUR" {print $1}' data/integrated_call_samples_v3.20130502.ALL.panel > data/EUR.txt
```

### Quick run

**File: 1000g.yaml**
```yaml
vcf_template: '${launchDir}/data/1000g.chr{CHR}'
keep: '${launchDir}/data/EUR.txt'
outdir: 'output'
ld_command: '--r-unphased ref-based cols=id,ref,alt'
prefix: 'EUR'
chroms: '20'
ld_window_kb: 1000
ld_window_r2: 0.01
min_col: 'UNPHASED_R'
ld_threads: 1
chunk_size_kb: 20000
overlap_size_kb: 1000
```

If `plink2` or `ldzip` are NOT available in your `$PATH`, export their paths:

```bash
export PLINK2=/path/to/plink2
export LDZIP=$(pwd)/LDZip/cpp/bin/ldzip
```

Run as follows:

```bash
nextflow run LDZip/pipelines/wholeGenomeLD/main.nf -params-file 1000g.yaml -resume
```

Go to Step 3 to test whether it worked correctly.

### Full run

To run on all chromosomes, update the YAML:

```yaml
chroms: '1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22'
ld_threads: 8
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

process.withName: ldPlink {
    cpus = params.ld_threads
}
executor.perCpuMemAllocation = true
```

Run as follows:

```bash
nextflow run LDZip/pipelines/wholeGenomeLD/main.nf -params-file 1000g.yaml -C slurm.config -resume
```

For other HPC environments, refer to the Nextflow executor [guidelines](https://www.nextflow.io/docs/latest/executor.html).

---

## Step 3: Verify Output

After successful completion, verify the LD matrix by querying two variants:

```r
library(LDZipMatrix)

# Load the LD matrix
ld <- LDZipMatrix("output/whole_genome/EUR")

# Query LD between two variants
fetchLD(ld, "20:64331475:C:T", "20:64333832:C:A")
# [1] 0.976378

# Benchmark query time
system.time(fetchLD(ld, "20:64331475:C:T", "20:64333832:C:A"))
#    user  system elapsed
#   0.003   0.000   0.010
```

---

## Parameters

| Parameter        | Description                                                                 | Example / Notes |
|------------------|-----------------------------------------------------------------------------|-----------------|
| `vcf_template`   | Template path to per-chromosome VCF files that you downloaded. `{CHR}` is replaced at runtime | `${launchDir}/data/1000g.chr{CHR}` |
| `keep`           | File with list of sample IDs to retain (generated above)              | EUR sample list |
| `outdir`         | Directory where LD outputs will be written                                | `output` |
| `ld_command`     | PLINK2 LD computation flags                                               | `--r-unphased ref-based cols=id,ref,alt` |
| `prefix`         | Prefix used for naming output files                                       | `EUR` |
| `chroms`         | Comma-separated list of chromosomes to process                            | `1–22` |
| `ld_window_kb`   | LD window size in kilobases                                               | `1000` (1 Mb) |
| `ld_window_r2`   | Minimum r² threshold for reporting LD pairs                               | `0.01` |
| `min_col`        | LD metric column to extract/store                                         | `UNPHASED_R` |
| `ld_threads`     | Number of threads used for LD computation                                 | `8` |
| `chunk_size_kb`  | *(Optional)* Chunk size in kb. Defaults to `2 × ld_window_kb`             | `2000` |
| `overlap_size_kb`| *(Optional)* Overlap size in kb. Defaults to `ld_window_kb`               | `1000` |
| `stage_chunk`    | *(Optional)* Stage intermediate chunk files to outdir (default: `false`)  | `true` |
| `stage_chr`      | *(Optional)* Stage per-chromosome files to outdir (default: `false`)      | `true` |