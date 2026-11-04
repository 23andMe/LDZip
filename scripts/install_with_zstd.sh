#!/usr/bin/env bash
set -euo pipefail

# Create installation directory
mkdir -p "$HOME/local"
cd "$HOME/local"

# Clone, build, and install zstd
git clone https://github.com/facebook/zstd.git
cd zstd
make -j"$(nproc)"
make install PREFIX="$HOME/local/zstd_install"

# Clone LDZip
cd "$HOME/local"
git clone https://github.com/23andMe/LDZip.git
cd LDZip/cpp

# Build LDZip against locally installed zstd
make build ZSTD_PREFIX="$HOME/local/zstd_install" -j"$(nproc)" V=1

# Use locally installed zstd at runtime
export LD_LIBRARY_PATH="$HOME/local/zstd_install/lib:${LD_LIBRARY_PATH:-}"

# Verify LDZip executable
ldzip --help

# Install R package dependencies
Rscript -e "p<-c('Rcpp','RSQLite','DBI'); p<-setdiff(p, rownames(installed.packages())); if(length(p)) install.packages(p, repos='https://cran.r-project.org')"

# Install R package
cd ../R
make