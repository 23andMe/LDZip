# Generating Whole-Genome LDZip Matrix from 1000 Genomes 

This tutorial walks through generating whole-genome linkage disequilibrium (LD) matrices from 1000 Genomes Project (1000G) data using the Nextflow pipeline provided in this repository. The same workflow can be applied to other large-scale datasets such as UK Biobank or TOPMed.

## Overview

The workflow consists of two main steps:

1. Download phased VCFs for each chromosome  
2. Run a Nextflow pipeline to compute a whole genome LDZip matrix  

---

## Prerequisites

- [`nextflow`](https://www.nextflow.io/docs/latest/install.html)
- [`plink2`](https://www.cog-genomics.org/plink/2.0/)
- Sufficient compute (multi-core recommended)
- Adequate storage for final/intermediate files

---

## Step 1: Download 1000G Phased VCFs

```bash
mkdir -p data

# download VCFs
for chr in {1..22}; do
  wget -O data/1000g.chr${chr}.vcf.gz ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased/CCDG_14151_B01_GRM_WGS_2020-08-05_chr${chr}.filtered.shapeit2-duohmm-phased.vcf.gz
  wget -O data/1000g.chr${chr}.vcf.gz.tbi ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased/CCDG_14151_B01_GRM_WGS_2020-08-05_chr${chr}.filtered.shapeit2-duohmm-phased.vcf.gz.tbi
done

# download sample panel (2504 samples, phase3)
wget https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel

```


---

## Step 2: Run Nextflow LD Pipeline for European Samples


- Clone the LDZip repository:

  ```bash
  git clone git@github.com:23andMe/LDZip.git
  ```

- Retrieve list of EUR samples

  ```bash
  awk '$3=="EUR" {print $1}' integrated_call_samples_v3.20130502.ALL.panel > data/EUR.txt
  ```

- Create YAML file of input parameters

  **File: 1000g.yaml**
  ```yaml
  vcf_template: '${launchDir}/data/1000g.chr{CHR}'
  keep: '${launchDir}/data/EUR.txt'
  outdir: 'output'
  ld_command: '--r-unphased ref-based cols=id,ref,alt'
  prefix: 'EUR'
  chroms: '20,21,22'
  ld_window_kb: 1
  ld_window_r2: 0.2
  min_col: 'UNPHASED_R'
  ld_threads: 1
  ```

- If `plink2` or `ldzip` are NOT available in your `$PATH`, export their paths before running the pipeline (might need to use `set -x` or `setenv` depending on the shell):

  ```bash
  export PLINK2=/path/to/plink2
  export LDZIP=$(pwd)/LDZip/cpp/bin/ldzip
  ```

- Run the pipeline:

  ```bash
  nextflow run LDZip/pipelines/wholeGenomeLD/main.nf -params-file 1000g.yaml -resume
  ```

  ---

- **Note**

  The above configuration runs locally for chromosomes 20-22 only (`chroms`), with a small LD window (`ld_window_kb: 1`) and a higher r² threshold (`ld_window_r2: 0.2`) so the pipeline runs quickly (typically ~5 minutes). This serves as a sanity check to ensure the inputs are correct and the workflow completes end-to-end with reasonable outputs.

  To run LDZip on all chromosomes, update these parameters to more reasonable values


  ```yaml
  chroms: '1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22'
  ld_window_kb: 1000
  ld_window_r2: 0.01
  ld_threads: 8
  ```

- **Running on HPC**

  If you are running with a higher resolution as above, you would want to use a HPC. E.g. if you were using a SLURM cluster you could pass the following config file to your command line

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

  ```bash
  nextflow run LDZip/pipelines/wholeGenomeLD/main.nf -params-file 1000g.yaml -C slurm.config -resume
  ```

  For other HPC environments, refer to the Nextflow executor [guidelines](https://www.nextflow.io/docs/latest/executor.html).
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