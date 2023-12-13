#!/usr/bin/env python3

import json
import logging
import anndata as ad
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


# Read in the configuration
with open("input_config.json", "r") as handle:
    config: dict = json.load(handle)
log("Read in configuration")
log(json.dumps(config, indent=4))

log("Reading")
adata = ad.read_h5ad("adata.h5ad")
log(str(adata))

# log-transform
log("Log transforming")
sc.pp.log1p(adata)

# Write out the normalized dataset
adata.write_h5ad("normalized.h5ad")

# Write out the updated configuration
log("Writing out configuration")
log(json.dumps(config, indent=4))
with open("output_config.json", "w") as handle:
    json.dump(config, handle, indent=4)
