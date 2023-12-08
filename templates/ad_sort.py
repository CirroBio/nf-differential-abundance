#!/usr/bin/env python3

import anndata as ad
from scipy.cluster import hierarchy

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
