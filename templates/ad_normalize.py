#!/usr/bin/env python3

import anndata as ad
import scanpy as sc

print("Reading")
adata = ad.read_h5ad("adata.h5ad")
print(adata)

# Normalize to proportions and log-transform
print("Normalizing")
sc.pp.normalize_per_cell(adata)
sc.pp.log1p(adata)

# Write out the normalized dataset
adata.write_h5ad("normalized.h5ad")
