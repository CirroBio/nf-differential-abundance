#!/usr/bin/env python3

import json
import anndata as ad
import scanpy as sc

# Read in the configuration
with open("input_config.json", "r") as handle:
    config: dict = json.load(handle)
print("Read in configuration")
print(json.dumps(config, indent=4))

print("Reading")
adata = ad.read_h5ad("adata.h5ad")
print(adata)

# log-transform
print("Log transforming")
sc.pp.log1p(adata)

# Write out the normalized dataset
adata.write_h5ad("normalized.h5ad")

# Write out the updated configuration
print("Writing out configuration")
print(json.dumps(config, indent=4))
with open("output_config.json", "w") as handle:
    json.dump(config, handle, indent=4)
