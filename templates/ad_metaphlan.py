#!/usr/bin/env python3

import json
import numpy as np
import pandas as pd
from anndata import AnnData
from pathlib import Path


# Scale values from -0.5 to 0.5
def scale_values(v: pd.Series):
    return -0.5 + (v - v.min()) / (v.max() - v.min())


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

# Make a list of the annotations provided
obs_sets = [
    dict(
        name=cname,
        path=f"obs/{cname}"
    )
    for cname in adata.obs.columns.values
]

# Filter by minimum abundance
min_abund = float("${params.min_abund}")
print(f"Minimum abundance threshold: {min_abund}")
mask = adata.to_df().max() > min_abund
print(f"Filtering to {mask.sum():,} / {mask.shape[0]:,} taxa")
adata = adata[:, mask]

# Start building the list of elements for visualization
config = dict()
for stats_fp in Path("corncob/").rglob("*.csv"):
    # Read the table
    df = pd.read_csv(stats_fp, index_col=0)
    # Make the order match
    df = df.reindex(index=adata.var_names)
    # Add the -log10(pvalue)
    df = df.assign(
        neg_log10_pvalue=-df["p_value"].apply(np.log10),
        sig_diff=df["p_value"] <= 0.05,
        mean_abund=adata.to_df().fillna(0).mean()
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
        volcano_varm = kw + "_volcano"
        ma_varm = kw + "_ma"
        adata.varm[volcano_varm] = (
            df
            .query("neg_log10_pvalue > 0.1")
            .reindex(
                columns=["est_coef", "neg_log10_pvalue"],
                index=adata.var_names
            )
            .apply(scale_values)
            .values
        )

        # Also make an MA plot
        adata.varm[ma_varm] = (
            df
            .query("neg_log10_pvalue > 0.1")
            .reindex(
                columns=["est_coef", "mean_abund"],
                index=adata.var_names
            )
            .apply(scale_values)
            .values
        )

        # Sort the annotations for this visualization
        obs_sets.sort(
            key=lambda d: d['name'] != kw
        )

        # Set up the annotations for the variables
        var_sets = [
            {
                "name": f"{kw}-Associated",
                "path": f"obsm/{name}/sig_diff"
            },
            {
                "name": "Mean Proportion (%)",
                "path": "obs/mean_proportion"
            }
        ]

        # Pick the top organism to select
        top_taxon = (
            df.apply(
                lambda r: r['mean_abund'] * r['neg_log10_pvalue'],
                axis=1
            )
            .sort_values(ascending=False)
            .index.values[0]
        )

        config[kw] = dict(
            title=f"Microbiome ~ {kw}",
            description=f"Summary of organisms associated with {kw}",
            obs_title=f"{tax_level.title()} Level Abundances",
            obs_sets=obs_sets,
            var_sets=var_sets,
            obs_type="sample",
            feature_type="species",
            volcano_varm=volcano_varm,
            volcano_title=f"{kw} - Volcano Plot",
            ma_varm=ma_varm,
            ma_title=f"{kw} - MA Plot",
            top_taxon=top_taxon
        )

# Write out the full data
adata.write_h5ad("metaphlan.h5ad", compression="gzip")

# Write out the list of variables for visualization
with open("config.json", "w") as handle:
    json.dump(config, handle, indent=4)
