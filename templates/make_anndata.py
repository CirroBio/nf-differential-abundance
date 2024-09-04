#!/usr/bin/env python3

import json
import numpy as np
import pandas as pd
from anndata import AnnData
from pathlib import Path
import logging
import plotly.express as px

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
    log(f"Preparing to add the {kw} table")
    log(df.head().to_csv())
    log(f"Using columns: {', '.join(cnames)}")

    # Make a scatterplot
    fig = px.scatter(
        data_frame=df.reset_index(),
        x=cnames[0],
        y=cnames[1],
        hover_name=df.index.name,
        hover_data=df.columns.values,
        title=kw
    )
    fig.update_layout(
        plot_bgcolor='rgba(0, 0, 0, 0)',
        paper_bgcolor='rgba(0, 0, 0, 0)'
    )
    fig.write_html(f"{kw}.html")

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

    adata.varm[kw] = df
    adata.uns[f"varm_cnames_{kw}"] = df.columns.values


# Read in global config elements
tax_level = "${params.tax_level}"

# Read in all of the locally staged files
dat = {
    fp.name.replace(".csv", ""): pd.read_csv(fp, index_col=0)
    for fp in Path(".").rglob("*.csv")
    # Ignore the corncob results for now
    if not fp.name.startswith(("mu.", "phi."))
}

# Find any samples which have 0 counts
zero_counts = (
    dat["counts"].index.values
    [dat["counts"].sum(axis=1) == 0]
)
if len(zero_counts) > 0:
    log("Dropping samples with 0 counts:")
    log("\\n".join(map(str, zero_counts)))
    dat["counts"] = dat["counts"].drop(index=zero_counts)
    dat["proportions"] = dat["proportions"].drop(index=zero_counts)

for kw, df in dat.items():
    log(kw)
    log(df.head().to_csv())
    log(f"Rows: {df.shape[0]:,}")
    log(f"Columns: {df.shape[1]:,}")

# For each of the statistical tests, read in the results
# and make any adjustments to column names which are needed
stats_name_map = {
    "PValue": "p_value",
    "pvalue": "p_value",
    "pval": "p_value",
    "we.ep": "p_value",
    "P.Value": "p_value",
    "pvalues": "p_value",
    "coef": "est_coef",
    "effect": "est_coef",
    "logFC": "est_coef",
    "lfc": "est_coef",
    "log2FoldChange": "est_coef"
}
stats = {
    stats_fp.name.replace(".csv", ""): (
        pd.read_csv(stats_fp, index_col=0)
        .rename(columns=stats_name_map)
    )
    for stats_fp in Path("stats/").rglob("*.csv")
}

for kw, df in stats.items():
    log(kw)
    log(df.head().to_csv())
    log(f"Rows: {df.shape[0]:,}")
    log(f"Columns: {df.shape[1]:,}")
    assert "p_value" in df.columns.values
    assert "est_coef" in df.columns.values

# Make a list of any variables which have stats results
all_vars = list(set([
    var_name
    for df in stats.values()
    for var_name in df.index
]))

# Make an AnnData object
adata = AnnData(
    (
        dat
        ["proportions"]
        .apply(lambda r: r / r.sum(), axis=1)
        .reindex(columns=all_vars)
    ),
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
        .apply(
            lambda cvals: (
                cvals
                if cvals.apply(type).nunique() == 1
                else cvals.fillna("").astype(str)
            )
        )
    ),
    var=(
        dat["taxonomy"]
        .assign(
            mean_proportion=dat["proportions"].mean(),
            mean_prevalence=(dat["proportions"] > 0).mean()
        )
        .reindex(index=all_vars)
    ),
    layers=dict(
        counts=(
            dat["counts"]
            .reindex(columns=all_vars)
        )
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

# Color taxa that exceed the fdr_cutoff
fdr_cutoff = float("${params.fdr_cutoff}")

# Start building the list of elements for visualization
config = dict()
for name, df in stats.items():
    # Name for the stat
    log(name)

    # Make the order match
    df = (
        df
        .reindex(index=adata.var_names)
    )

    # Drop any columns which are all NaN
    df = df.dropna(axis=1, how="all")

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

    # Save the whole table
    adata.varm[name] = df
    adata.uns[f"varm_cnames_{name}"] = df.columns.values

    # For the mu., make a volcano plot
    if name.startswith("mu.") or "${params.method}" != "corncob":

        # The statistical method used will impact the:
        # - parsing of metric name from file name
        # - variable used to show effect size
        effect_cname = "est_coef"
        effect_title = "Fold Change (log10)"
        if "${params.method}" == "corncob":
            kw = name[3:]
        else:
            kw = name

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
                "name": f"{kw}-Associated",
                "path": f"obsm/{name}/sig_diff"
            },
            {
                "name": "p-value",
                "path": f"obsm/{name}/p_value"
            },
            {
                "name": effect_title,
                "path": f"obsm/{name}/{effect_cname}"
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
