#!/usr/bin/env python3

import json
import pandas as pd

levels = [
    "kingdom",
    "phylum",
    "class",
    "order",
    "family",
    "genus",
    "species"
]


def read_table(fp, **kwargs) -> pd.DataFrame:
    print(f"Reading in: {fp}")
    df = pd.read_csv(fp, **kwargs)
    print(df.to_csv())
    return df


if __name__ == "__main__":

    # Read in the table of abundances
    print("Parsing params.read_kwargs")
    read_kwargs = json.loads('${params.read_kwargs}')
    counts = read_table("${abund}", **read_kwargs)
    counts.index.names = ['sample']

    # Read in the table of metadata
    metadata = read_table("metadata.csv", index_col=0)
    metadata.index.names = ['sample']

    # Only write out samples which are shared
    shared_samples = list(
        set(counts.index.values) &
        set(metadata.index.values)
    )
    shared_samples.sort()
    print(f"Samples with metadata: {metadata.shape[0]:,}")
    print(f"Samples with counts: {counts.shape[0]:,}")
    print(f"Samples shared: {len(shared_samples):,}")

    assert len(shared_samples) > 1

    # Set up a dummy taxonomy
    taxonomy = pd.DataFrame([
        {
            "name": feature,
            **{
                level.title(): feature
                for level in levels
            }
        }
        for feature in counts.columns.values
    ]).set_index("name")
    print("Dummy taxonomy")
    print(taxonomy.to_csv())

    # Write out files
    metadata.reindex(index=shared_samples).to_csv("samplesheet.csv")
    counts = counts.reindex(index=shared_samples)
    counts.to_csv("counts.csv")
    counts.apply(lambda r: r / r.sum(), axis=1).to_csv("proportions.csv")
    taxonomy.to_csv("taxonomy.csv")
