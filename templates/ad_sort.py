#!/usr/bin/env python3

import json
import anndata as ad
from scipy.cluster import hierarchy

# Read in the configuration
with open("input_config.json", "r") as handle:
    config: dict = json.load(handle)
print("Read in configuration")
print(json.dumps(config, indent=4))

print("Reading")
adata = ad.read_h5ad("adata.h5ad")
print(adata)


# Sort features and observations
def sort(adata: ad.AnnData):
    return adata[
        adata.obs_names[
            hierarchy.leaves_list(
                hierarchy.linkage(
                    adata.X,
                    method="ward"
                )
            )
        ]
    ].copy()


print("Sorting observations")
adata = sort(adata)
print("Sorting features")
adata = sort(adata.copy().T).T

# Write out the sorted dataset
adata.write_h5ad("sorted.h5ad")

# Write out the updated configuration
print("Writing out configuration")
print(json.dumps(config, indent=4))
with open("output_config.json", "w") as handle:
    json.dump(config, handle, indent=4)
