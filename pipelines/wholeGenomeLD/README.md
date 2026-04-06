# Nextflow Whole Genome LD Pipeline

This pipeline manages the scatter gather to generate a compressed whole genome LD matrix from PLINK2 genotype data. It uses PLINK2 to calculate the LD for each chromosome, uses `ldzip` to compress each matrix, block-diagonal concatenates them to generate the whole genome matrix, and finally builds and SQLITE database for fast querying with variant IDs.

**For detailed parameter descriptions and validation rules, run:**
```bash
nextflow run main.nf --help
```

---

## Run Quick Test

| Note: Make sure your C++ binary `LDZip/cpp/bin/ldzip` is available in `$PATH` for the Nextflow pipelines to work. Otherwise, provide the path to the binary as `env.LDZIP` in `nextflow.config`

```bash
./main.nf -params-file tests/yaml/unit1.yaml 
```
---

## Usage Parameters

### Input parameters

**Genotype Input** (provide ONE of the following):
- `--pfile_template` : Template for PLINK2 `.pfile` inputs (use `{CHR}` as placeholder, e.g. `/home/data/g1k.chr{CHR}`).
- `--vcf_template` : Template for VCF files (use `{CHR}` as placeholder, e.g. `/home/data/g1k.chr{CHR}`). VCF files will be automatically converted to pgen format.

**Other Input Options**:
- `--chroms` : Comma-separated list of chromosomes to process (e.g. `1,2,X,Y`). Defaults to autosomes,X,Y and MT.
- `--extract` : File with variant IDs to include.
- `--exclude` : File with variant IDs to exclude.
- `--keep` : File with sample IDs to subset.

### PLINK2 LD parameters
- `--ld_command` : Command to pass to PLINK2 (default: `--r-phased ref-based cols=id,ref,alt,dprime`).
- `--ld_threads` : Number of CPU threads for LD calculations (default: `4`).
- `--ld_window_kb` : Window size in kb for LD computation (default: `1000`).
- `--ld_window_r2` : Minimum r² threshold for reporting LD pairs (default: `0.01`).

### Chunking parameters
- `--chunk_size_kb` : Size of each chunk in kilobases. If not specified, defaults to `2 × ld_window_kb`.
- `--overlap_size_kb` : Overlap between adjacent chunks in kilobases. If not specified, defaults to `ld_window_kb`.

### LD Compression parameters
- `--bits` : Number of bits for LD quantization (`8`, `16`, `32`, or `99` for raw float).
- `--min` : Minimum absolute LD value to retain (default: `0.0`).
- `--min_col` : Column to apply the above `min` filter on (default: `PHASED_R`).

### Output parameters
- `--outdir` : Output directory (default: current directory).
- `--prefix` : Prefix for final output files (default: `out`).
- `--stage_pgen` : Stage converted pgen files to outdir (only when using `--vcf_template`, default: `false`).
- `--stage_plink` : Stage intermediate PLINK LD files to outdir (default: `false`).
- `--stage_chunk` : Stage compressed LDZip chunk files to outdir (default: `false`).
- `--stage_chr` : Stage per-chromosome concatenated LDZip files to outdir (default: `false`).

---
