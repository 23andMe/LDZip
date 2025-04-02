# LDZip C++ Library

- The source code can be compiled to produce a standalone `ldzip` binary for efficient compression of PLINK LD matrices.  
- In addition, the library provides C++ headers exposing the `LDZip` class, which can be directly included and used by other C++ projects.  
- This same C++ codebase has been leveraged to build the `LDZipMatrix` R API, demonstrating how the core compression and random-access functionality can be reused across different environments.


## Install `ldzip` binary

The C++ binary is ONLY necessary if you are compressing a plink LD matrix. If you need to access an existing compressed matrix, it should be sufficient to install the R package `LDZipMatrix` (see below)

```bash
make build
```

## Run Quick Test 

```bash
make quick-test
```

## Usage

Once built, the `ldzip` binary is located at `cpp/bin/ldzip`.

It supports the following main operations, each with its own subcommands and options:

1. **compress** – Convert PLINK LD matrices (binary or tabular) into the efficient `.ldzip` format.  
2. **decompress** – Restore a `.ldzip` archive back into PLINK tabular or binary LD matrices.  
3. **concat** – Merge multiple `.ldzip` chunks into a single combined archive.  
4. **filter** – Filter a `.ldzip` file based on 0-based indices to another `.ldzip` file.  

---

### Compress


The `compress` operation has two subcommands, depending on the input format:

- **`plinkSquare`** – Use this when the input is a **square PLINK binary LD matrix** (usually with a `.unphased.vcor1.bin` suffix) file plus a SNP list).  
- **`plinkTabular`** – Use this when the input is a **PLINK tabular LD file** (text-based file usually with a `.vcor` suffix).  file plus a SNP list).  

Both subcommands share the same options:

| Option          | Description                                                                                                                               | Required / Default |
|-----------------|-------------------------------------------------------------------------------------------------------------------------------------------|--------------------|
| `--ld_file`     | Path to the input LD matrix file. For `plinkSquare`, this is the binary `.ld` file; for `plinkTabular`, this is the tabular LD file.       | **Required**       |
| `--snp_file`    | Path to the SNP list file corresponding to the LD matrix. This ensures correct indexing of variants.                                       | **Required**       |
| `--output_prefix` | Prefix for the compressed `.ldzip` output files. The compressor generates multiple files with this prefix (e.g. `.x.bin`, `.i.bin`, `io.bin`, `.p.bin`, `.meta.json`). | **Required**       |
| `--bits`        | Compression precision level. Supports `8`, `16`, `32`, or `99` (no quantization). Lower bit-widths reduce file size at the cost of precision. | Default: *8*       |
| `--min`         | Minimum absolute LD threshold. Pairs with LD values below this cutoff are discarded to save space.                                       | Default: *1e-4*    |
| `--format`      | Matrix storage format. Can be `upper` (store only the upper triangle, halving space) or `full` (store the entire matrix).                  | Default: *full*    |



**Example Usage: Plink Binary Input**

```
bin/ldzip compress plinkSquare  \
    --ld_file ../assets/unit.bin \
    --snp_file ../assets/unit.bin.vars \
    --output_prefix test \
    --bits 8 \
    --min 1e-4 
```

---


### Decompress


The `decompress` operation has the following options:

| Option          | Description                                                                                                                               | Required / Default |
|-----------------|-------------------------------------------------------------------------------------------------------------------------------------------|--------------------|
| `--input_prefix`     | Path to the input compressed LD matrix file (expects `.x.bin`, `.i.bin`, `io.bin`,  `.p.bin`, `.meta.json`)      | **Required**       |
| `--output_prefix` | Prefix for the uncompressed output files | **Required**       |
| `--type`        | Output Type (tabular or binary)  | Default: *binary*       |



**Example Usage: Plink Binary Input**

```
bin/ldzip decompress \
    --input_prefix test \
    --type binary \
    --output_prefix test_decompressed 
```


---


### Filter


The `filter` operation has the following options:

| Option                | Description                                                                                                                               | Required / Default |
|-----------------------|-------------------------------------------------------------------------------------------------------------------------------------------|--------------------|
| `--input_prefix`      | Path to the input compressed LD matrix file (expects `.x.bin`, `.i.bin`, `io.bin`, `.p.bin`, `.meta.json`)      | **Required**       |
| `--output_prefix`     | Prefix for the filtered compressed LD matrix file | **Required**       |
| `--range`             | Range of indices to filter (0-based START-END)    |       |
| `--keep`              | File with list of indices (0-based)               | Default: *binary*       |



**Example Usage: Plink Binary Input**

```
bin/ldzip filter \
    --input_prefix test \
    --range 5-15 \
    --output_prefix test_filtered 
```
