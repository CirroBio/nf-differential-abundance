#!/bin/bash
set -e

Rscript ${script} counts.tsv metadata.tsv ${params.method}.results.csv 2>&1 \
    | tee run_${params.method}.log

# Convert the TSV output to CSV
sed -i 's/\\t/,/g' ${params.method}.results.csv
