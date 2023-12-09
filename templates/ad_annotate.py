#!/usr/bin/env python3

import json
import anndata as ad
import scanpy as sc

adata = ad.read_h5ad("adata.h5ad")
print(adata)

# Read in the configuration
with open("config.json", "r") as handle:
    config: dict = json.load(handle)
print("Read in configuration")
print(json.dumps(config, indent=4))


def update_ordination_path(kw, title):
    path = f"obsm/X_{kw}"
    print(f"Setting ordination to use {path}")
    for kw in config.keys():
        config[kw]["ordination_path"] = path
        config[kw]["ordination_title"] = title


# Run t-SNE, UMAP, and PCA ordination
# Use a try/except pattern in case individual ones fail
print("Ordinating with PCA")
sc.tl.pca(adata)
update_ordination_path("pca", "PCA")

#  t-SNE
print("Ordinating with t-SNE")
metric = 'braycurtis'
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


# UMAP
print("Ordinating with UMAP")
try:
    sc.tl.umap(adata)
    update_ordination_path("umap", "UMAP")
except Exception as e:
    print("Skipping UMAP")
    print(str(e))

# Write out the annotated dataset
adata.write_h5ad("annotated.h5ad")

# Write out the updated configuration
print("Writing out configuration")
print(json.dumps(config, indent=4))
with open("config.json", "w") as handle:
    json.dump(config, handle, indent=4)
