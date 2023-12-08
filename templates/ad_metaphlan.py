#!/usr/bin/env python3

import pandas as pd
from anndata import AnnData
from pathlib import Path

# Read in all of the locally staged files
dat = {
    fp.name.replace(".csv", ""): pd.read_csv(fp, index_col=0)
    for fp in Path(".").rglob("*.csv")
    # Ignore the corncob results for now
    if not fp.name.startswith(("mu.", "phi."))
}

assert "mu.CRC" not in dat

for kw, df in dat.items():
    print(kw)
    print(df.head())
    print(df.shape)

# Make an AnnData object
ad = AnnData(
    dat["proportions"],
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
    var=dat["taxonomy"],
    layers=dict(
        counts=dat["counts"]
    )
)

# Add the stats
variables = []
for stats_fp in Path("corncob/").rglob("*.csv"):
    df = pd.read_csv(stats_fp, index_col=0)
    name = stats_fp.name.replace(".csv", "")
    print(name)
    ad.varm[name] = df.reindex(index=ad.var_names)
    print(ad.varm[name])
    variables.append(name)

# Write out the full data
ad.write_h5ad("metaphlan.h5ad", compression="gzip")

# Write out the list of variables for visuzliation
with open("variables.txt", "w") as handle:
    handle.write("\n".join(variables))
