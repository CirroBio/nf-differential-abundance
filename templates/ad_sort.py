#!/usr/bin/env python3

import json
import anndata as ad
from scipy.cluster import hierarchy
import logging

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


# Read in the configuration
with open("input_config.json", "r") as handle:
    config: dict = json.load(handle)
log("Read in configuration")
log(json.dumps(config, indent=4))

log("Reading")
adata = ad.read_h5ad("adata.h5ad")
log(str(adata))


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


log("Sorting observations")
adata = sort(adata)
log("Sorting features")
adata = sort(adata.copy().T).T

# Write out the sorted dataset
adata.write_h5ad("sorted.h5ad")

# Write out the updated configuration
log("Writing out configuration")
log(json.dumps(config, indent=4))
with open("output_config.json", "w") as handle:
    json.dump(config, handle, indent=4)
