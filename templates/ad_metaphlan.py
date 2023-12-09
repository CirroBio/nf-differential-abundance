#!/usr/bin/env python3

import json
import numpy as np
import pandas as pd
from anndata import AnnData
from pathlib import Path

# Read in global config elements
tax_level = "${params.tax_level}"

# Read in all of the locally staged files
dat = {
    fp.name.replace(".csv", ""): pd.read_csv(fp, index_col=0)
    for fp in Path(".").rglob("*.csv")
    # Ignore the corncob results for now
    if not fp.name.startswith(("mu.", "phi."))
}

for kw, df in dat.items():
    print(kw)
    print(df.head())
    print(df.shape)

# Make an AnnData object
adata = AnnData(
    dat["proportions"].apply(lambda r: r / r.sum(), axis=1),
    obs=(
        # Make sure that we're indexing on the sample column
        dat["samplesheet"]
        .reset_index()
        .set_index("sample")
        # Only take the samples which we have data for
        .reindex(
            index=dat["proportions"].index
        )
        .drop(
            columns=["file"]
        )
    ),
    var=(
        dat["taxonomy"]
        .assign(
            mean_proportion=dat["proportions"].mean()
        )
    ),
    layers=dict(
        counts=dat["counts"]
    )
)

# Filter by minimum abundance
min_abund = float("${params.min_abund}")
print(f"Minimum abundance threshold: {min_abund}")
mask = adata.to_df().max() > min_abund
print(f"Filtering to {mask.sum():,} / {mask.shape[0]:,} taxa")
adata = adata[:, mask]

# Start building the list of elements for visualization
config = []
for stats_fp in Path("corncob/").rglob("*.csv"):
    # Read the table
    df = pd.read_csv(stats_fp, index_col=0)
    # Make the order match
    df = df.reindex(index=adata.var_names)
    # Add the -log10(pvalue)
    df = df.assign(
        neg_log10_pvalue=-df["p_value"].apply(np.log10),
        sig_diff=df["p_value"] <= 0.05
    )
    # Name for the stat
    name = stats_fp.name.replace(".csv", "")
    print(name)
    # Save the whole table
    adata.varm[name] = df
    print(adata.varm[name])

    # For the mu., make a volcano plot
    if name.startswith("mu."):
        kw = name[3:]
        varm_volcano = kw + "_volcano"
        varm_ma = kw + "_ma"
        adata.varm[varm_volcano] = (
            df
            .query("neg_log10_pvalue > 0.1")
            .reindex(
                columns=["est_coef", "neg_log10_pvalue"],
                index=adata.var_names
            )
            .values
        )

        # Also make an MA plot
        adata.varm[varm_ma] = (
            df
            .query("neg_log10_pvalue > 0.1")
            .reindex(
                columns=["est_coef"],
                index=adata.var_names
            )
            .assign(
                mean_abund=adata.to_df().fillna(0).mean()
            )
            .values
        )

        config[kw] = dict(
            title=f"Microbiome ~ {kw}",
            description=f"Summary of organisms associated with {kw}",
            obs_title=f"{tax_level.title()} Level Abundances (observations)"
        )

# Write out the full data
adata.write_h5ad("metaphlan.h5ad", compression="gzip")

# Write out the list of variables for visualization
with open("config.json", "w") as handle:
    json.dump(config, handle, indent=4)
