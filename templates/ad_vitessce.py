#!/usr/bin/env python3

import anndata as ad

adata = ad.read_h5ad("adata.h5ad")
# Write out in both orientations
adata.write_zarr("obs.zarr")
adata.T.write_zarr("var.zarr")
