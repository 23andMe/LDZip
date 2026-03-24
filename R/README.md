# LDZipMatrix R package

This R package provides convenient access to compressed LD matrices, enabling efficient random access to LD values and related operations. Below are brief descriptions of some core functions. Please refer to the R help pages (`?function_name`) for detailed usage and arguments. Ensure that `roxygen2` is installed for documentation and `NAMESPACE` generation.



## Install R Package

```
git clone git@github.com:23andMe/LDZip.git
cd LDZip/R
make install
```

## Run Quick Test
| **Note**: You must have `tinytest` installed for the `quick-test` to work
```
make quick-test
```

## Usage

Please see the [R manual](man/LDZipMatrix.pdf) for more details and updated documentation. Below is a summary of the most used functions provided in the library

### - `LDZipMatrix()`

This function constructs an external pointer to a C++ `LDZipMatrix` object from a compressed file on disk, with which one can perform random access on the LD data as shown following this snippet.

```
library(LDZipMatrix)
ld = LDZipMatrix(file.path(system.file("extdata", package = "LDZipMatrix"), "g1k.chr22.ldzip"))
```
---
### - `buildIndex()`
Builds an SQLITE index file for querying variants by position or variant IDs.
- This step is necessary for using `fetchLD`/`getNeighbors` with variantIDs. 
- Integer indices will always work without building the SQLITE index
```
library(LDZipMatrix)
ld = LDZipMatrix(file.path(system.file("extdata", package = "LDZipMatrix"), "g1k.chr22.ldzip"))
buildIndex(ld)
```
---
### - `fetchLD(rows, columns)`
Retrieve entries from a compressed LD matrix by specifying rows and columns.  
- Supports queries by variant indices or variant IDs.  
- See `?fetchLD` for detailed usage and examples.
```
library(LDZipMatrix)
ld = LDZipMatrix(file.path(system.file("extdata", package = "LDZipMatrix"), "g1k.chr22.ldzip"))
fetchLD(ld, 2, 191, type = c("PHASED_R", "DPRIME"))
buildIndex(ld)
fetchLD(ld, "rs587755077", "rs587631919", type = "PHASED_R")
```
---

### - `getNeighbors(variant, type, abs_threshold)`
Find neighbors of a given variant within the LD matrix above a specified threshold for a given statistic type.  
- Returns variants with high LD greater than the cutoff and within a boundary
- See `?getNeighbors` for detailed usage and examples.
```
library(LDZipMatrix)
ld = LDZipMatrix(file.path(system.file("extdata", package = "LDZipMatrix"), "g1k.chr22.ldzip"))
buildIndex(ld) 
getNeighbors(ld, "rs587755077", type="PHASED_R", abs_threshold = sqrt(0.8), genomic_length=500000)
```

