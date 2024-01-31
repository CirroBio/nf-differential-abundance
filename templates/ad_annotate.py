#!/usr/bin/env python3

import anndata as ad
import json
import logging
import pandas as pd
import scanpy as sc

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        logging.FileHandler("${task.process}.log"),
        logging.StreamHandler()
    ]
)


def log(s: str):
    for line in s.split("\\n"):
        logging.info(line)


adata = ad.read_h5ad("adata.h5ad")
log(str(adata))

# Read in the configuration
with open("input_config.json", "r") as handle:
    config: dict = json.load(handle)
log("Read in configuration")
log(json.dumps(config, indent=4))


# Scale values from -0.5 to 0.5
def scale_values(v: pd.Series):
    # Do not divide by zero for invariant data
    span = v.max() - v.min()
    return -0.5 + (v - v.min()) / (span if span > 0 else 1)


def update_ordination_path(kw, title):
    path = f"obsm/X_{kw}"
    # Skip if there aren't enough dimensions to plot
    if adata.obsm[f"X_{kw}"].shape[1] < 2:
        log(f"Skipping {kw} - not enough dimensions")

    # Scale the values
    adata.obsm[f"X_{kw}"] = (
        pd.DataFrame(adata.obsm[f"X_{kw}"])
        .apply(scale_values)
        .values
    )

    log(f"Setting ordination to use {path}")
    log(str(adata.obsm[f"X_{kw}"]))
    for kw in config.keys():
        config[kw]["ordination_path"] = path
        config[kw]["ordination_title"] = title


# Run t-SNE, UMAP, and PCA ordination
# Use a try/except pattern in case individual ones fail
log("Ordinating with PCA")
sc.tl.pca(adata)
update_ordination_path("pca", "PCA")

# UMAP
log("Ordinating with UMAP")
metric = 'braycurtis'
try:
    sc.pp.neighbors(adata, use_rep='X', metric=metric)
    sc.tl.umap(adata)
    update_ordination_path("umap", "UMAP")
except Exception as e:
    log("Skipping UMAP")
    log(str(e))

#  t-SNE
log("Ordinating with t-SNE")
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
    log("Skipping t-SNE")
    log(str(e))

# Write out the annotated dataset
adata.write_h5ad("annotated.h5ad")

# Write out CSVs with each element
adata.obs.to_csv("output.obs.csv")
adata.var.to_csv("output.var.csv")
for kw in adata.varm:
    pd.DataFrame(
        adata.varm[kw],
        index=adata.var_names,
        columns=adata.uns.get(
            f"varm_cnames_{kw}",
            range(adata.varm[kw].shape[1])
        )
    ).to_csv(f"output.{kw}.csv")
for layer in adata.layers:
    adata.to_df(layer).to_csv(f"output.{layer}.csv")

# Write out the updated configuration
log("Writing out configuration")
log(json.dumps(config, indent=4))
with open("output_config.json", "w") as handle:
    json.dump(config, handle, indent=4)
