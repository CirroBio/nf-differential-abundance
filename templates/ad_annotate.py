#!/usr/bin/env python3

import anndata as ad
import scanpy as sc
from scipy.cluster import hierarchy

adata = ad.read_h5ad("adata.h5ad")
print(adata)

# Run t-SNE, UMAP, and PCA ordination
print("Ordinating with t-SNE")
metric = 'braycurtis'
sc.pp.neighbors(adata, use_rep='X', metric=metric)
sc.tl.tsne(
    adata,
    use_rep='X',
    perplexity=min(adata.n_obs - 1, 30),
    metric=metric
)
print("Ordinating with UMAP")
sc.tl.umap(adata)
print("Ordinating with PCA")
sc.tl.pca(adata)

# Write out the annotated dataset
adata.write_h5ad("annotated.h5ad")
