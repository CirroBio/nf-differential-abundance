#!/usr/bin/env python3

import json
import anndata as ad
import scanpy as sc
import pandas as pd

adata = ad.read_h5ad("adata.h5ad")
print(adata)

# Read in the configuration
with open("input_config.json", "r") as handle:
    config: dict = json.load(handle)
print("Read in configuration")
print(json.dumps(config, indent=4))


# Scale values from -0.5 to 0.5
def scale_values(v: pd.Series):
    # Do not divide by zero for invariant data
    span = v.max() - v.min()
    return -0.5 + (v - v.min()) / (span if span > 0 else 1)


def update_ordination_path(kw, title):
    path = f"obsm/X_{kw}"
    # Skip if there aren't enough dimensions to plot
    if adata.obsm[f"X_{kw}"].shape[1] < 2:
        print(f"Skipping {kw} - not enough dimensions")

    # Scale the values
    adata.obsm[f"X_{kw}"] = (
        pd.DataFrame(adata.obsm[f"X_{kw}"])
        .apply(scale_values)
        .values
    )

    print(f"Setting ordination to use {path}")
    print(adata.obsm[f"X_{kw}"])
    for kw in config.keys():
        config[kw]["ordination_path"] = path
        config[kw]["ordination_title"] = title


# Run t-SNE, UMAP, and PCA ordination
# Use a try/except pattern in case individual ones fail
print("Ordinating with PCA")
sc.tl.pca(adata)
update_ordination_path("pca", "PCA")

# UMAP
print("Ordinating with UMAP")
metric = 'braycurtis'
try:
    sc.pp.neighbors(adata, use_rep='X', metric=metric)
    sc.tl.umap(adata)
    update_ordination_path("umap", "UMAP")
except Exception as e:
    print("Skipping UMAP")
    print(str(e))

#  t-SNE
print("Ordinating with t-SNE")
try:
    sc.pp.neighbors(adata, use_rep='X', metric=metric)
    sc.tl.tsne(
        adata,
        use_rep='X',
        perplexity=min(adata.n_obs - 1, 30),
        metric=metric
    )
    update_ordination_path("tsne", "t-SNE")
except Exception as e:
    print("Skipping t-SNE")
    print(str(e))

# Write out the annotated dataset
adata.write_h5ad("annotated.h5ad")

# Write out the updated configuration
print("Writing out configuration")
print(json.dumps(config, indent=4))
with open("output_config.json", "w") as handle:
    json.dump(config, handle, indent=4)
