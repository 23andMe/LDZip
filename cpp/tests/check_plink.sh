#!/bin/bash
if [ -z "$PLINK2" ]; then
    echo "Please export PLINK2=/path/to/plink2" >&2
    exit 1
fi