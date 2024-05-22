#!/usr/bin/env python3

import pandas as pd

# Read in the counts table
counts = pd.read_csv("inputs/counts.csv", index_col=0)

# Read in the metadata table
metadata = pd.read_csv("inputs/metadata.csv", index_col=0)

# Only keep the single column of the metadata which contains the parameter of interest
category = "${params.formula}"
assert category in metadata.columns, f"'{category}' not found in metadata columns"

metadata = metadata.reindex(columns=[category])

# Drop any missing values
metadata = metadata.dropna()
print(f"There are {metadata.shape[0]:,} samples with information for '{category}'")
assert metadata.shape[0] > 0, f"No samples with information for '{category}'"

# Align the order of rows with the counts
counts = counts.loc[metadata.index]

# Write out as TSV
counts.T.to_csv("counts.tsv", sep="\t")
metadata.to_csv("metadata.tsv", sep="\t")
