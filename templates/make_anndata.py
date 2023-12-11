#!/usr/bin/env python3

import json
import numpy as np
import pandas as pd
from anndata import AnnData
from pathlib import Path


# Scale values from -0.5 to 0.5
def scale_values(v: pd.Series):
    # Do not divide by zero for invariant data
    span = v.max() - v.min()
    return -0.5 + (v - v.min()) / (span if span > 0 else 1)


# Add a table to the AnnData object, with scaled values
def add_varm(
    adata: AnnData,
    kw: str,
    df: pd.DataFrame,
    cnames: list
):
    print(f"Preparing to add the {kw} table")
    print(df.to_csv())
    print(f"Using columns: {', '.join(cnames)}")

    # Filter on the 'show' column
    assert "show" in df.columns.values
    df = df.query("show")
    assert df.shape[0] > 0
    assert df.dropna().shape[0] > 0, df.to_csv()

    # Just get the columns of interest
    df = df.reindex(
        columns=cnames,
        index=adata.var_names
    )
    assert df.shape[0] > 0
    assert df.dropna().shape[0] > 0, df.to_csv()

    # Scale the values
    df = df.apply(scale_values)
    assert df.shape[0] > 0
    assert df.dropna().shape[0] > 0, df.to_csv()

    adata.varm[kw] = df.values


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
        .rename(columns=dict(sample_id="sample"))
        .set_index("sample")
        # Only take the samples which we have data for
        .reindex(
            index=dat["proportions"].index
        )
        .drop(
            columns=(
                ["file"]
                if "file" in dat["samplesheet"].columns.values
                else []
            )
        )
    ),
    var=(
        dat["taxonomy"]
        .assign(
            mean_proportion=dat["proportions"].mean(),
            mean_prevalence=(dat["proportions"] > 0).mean()
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

# Color taxa that exceed the fdr_cutoff
fdr_cutoff = float("${params.fdr_cutoff}")


# Start building the list of elements for visualization
config = dict()
for stats_fp in Path("stats/").rglob("*.csv"):

    # Read the table
    df = pd.read_csv(stats_fp, index_col=0)

    # Make the order match
    df = df.reindex(index=adata.var_names)

    # Add the -log10(pvalue)
    df = df.assign(
        neg_log10_pvalue=-df["p_value"].apply(np.log10),
        sig_diff=df["p_value"] <= fdr_cutoff,
        mean_abund=adata.to_df().fillna(0).mean()
    )

    # Only show features in the volcano and MA plots which
    # have a neg_log10_pvalue > 0.1,
    # unless there are < 10 features under the threshold,
    # in which case all should be shown
    df = df.assign(
        show=(
            df["neg_log10_pvalue"] > 0.1
            if (df["neg_log10_pvalue"] > 0.1).sum() > 10
            else True
        )
    )

    # Name for the stat
    name = stats_fp.name.replace(".csv", "")
    print(name)
    # Save the whole table
    adata.varm[name] = df
    print(adata.varm[name])

    # For the mu., make a volcano plot
    if name.startswith("mu.") or "${params.method}" == "wilcoxon":

        # The statistical method used will impact the:
        # - parsing of metric name from file name
        # - variable used to show effect size
        if "${params.method}" == "wilcoxon":
            kw = name
            effect_cname = "lfc"
            effect_title = "Fold Change (log10)"
        else:
            kw = name[3:]
            effect_cname = "est_coef"
            effect_title = "Estimated Coefficient of Association"

        volcano_varm = kw + "_volcano"
        add_varm(adata, volcano_varm, df, [effect_cname, "neg_log10_pvalue"])

        # Also make an MA plot
        ma_varm = kw + "_ma"
        add_varm(adata, ma_varm, df, [effect_cname, "mean_abund"])

        # Sort the annotations for this visualization
        obs_sets.sort(
            key=lambda d: f"{(d['name'] != kw) and (kw.startswith(d['name']) is False)}-{d['name']}"
        )

        # Set up the annotations for the variables
        # (note: obsm is used instead of varm because
        # the data will ultimately be transposed for viewing)
        var_sets = [
            {
                "name": "p-value",
                "path": f"obsm/{name}/p_value"
            },
            {
                "name": effect_title,
                "path": f"obsm/{name}/{effect_cname}"
            },
            {
                "name": f"{kw}-Associated",
                "path": f"obsm/{name}/sig_diff"
            },
            {
                "name": "Average Abundance (%)",
                "path": "obs/mean_proportion"
            },
            {
                "name": "Proportion Detected (%)",
                "path": "obs/mean_prevalence"
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
adata.write_h5ad("${params.tool}.h5ad", compression="gzip")

# Write out the list of variables for visualization
with open("output_config.json", "w") as handle:
    json.dump(config, handle, indent=4)

# Write out a short summary
with open("${params.tool}.h5ad.txt", "w") as handle:
    handle.write(str(adata))
