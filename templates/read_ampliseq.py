#!/usr/bin/env python3

import pandas as pd
from typing import Tuple
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


mapping = dict(
    d="kingdom",
    p="phylum",
    c="class",
    o="order",
    f="family",
    g="genus",
    s="species"
)


def read_table(fp, **kwargs) -> pd.DataFrame:
    log(f"Reading in: {fp}")
    df = pd.read_csv(fp, **kwargs)
    log(df.head().to_csv())
    return df


def read_taxonomy() -> pd.DataFrame:
    taxonomy = read_table(
        "taxonomy.tsv",
        sep="\t",
        index_col=0
    )
    msg = "Expected to see 'Taxon' in columns"
    assert "Taxon" in taxonomy.columns.values, msg

    # Reformat the taxonomy
    taxonomy = pd.DataFrame([
        dict(
            asv=asv,
            rank=rank,
            name=name
        )
        for asv, tax_string in taxonomy["Taxon"].items()
        for rank, name in parse_tax_string(tax_string)
    ]).pivot(
        index='asv',
        columns='rank',
        values='name'
    )
    log("Reformatted taxonomy")
    log(taxonomy.head().to_csv())
    return taxonomy


def parse_tax_string(tax_string: str):

    for taxon in tax_string.split(";"):
        taxon = taxon.strip()
        assert taxon[1:3] == "__", taxon
        assert taxon[0] in mapping, taxon
        yield mapping[taxon[0]], taxon[3:]


def read_metadata() -> pd.DataFrame:

    return (
        read_table("metadata.csv")
        .set_index("sample")
    )


def read_counts() -> pd.DataFrame:

    return (
        read_table(
            "feature-table.tsv",
            sep="\t",
            skiprows=1,
            index_col=0
        )
        .applymap(int)
    )


def merge_counts(
    counts: pd.DataFrame,
    taxonomy: pd.DataFrame
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    tax_level = "${params.tax_level}"
    log(f"Merging counts by {tax_level}")
    assert tax_level in taxonomy.columns.values, taxonomy.columns.values

    counts = (
        counts
        .reindex(
            index=taxonomy[tax_level].dropna().index
        )
        .groupby(
            taxonomy[tax_level]
        )
        .sum()
        .sort_index()
    )
    log(counts.head().to_csv())

    # Do the same thing for the taxonomy
    taxonomy = (
        taxonomy
        .assign(name=taxonomy[tax_level])
        .reindex(index=taxonomy[tax_level].dropna().index)
        .groupby('name')
        .head(1)
        .set_index('name')
        .sort_index()
    )
    return counts, taxonomy


if __name__ == "__main__":

    taxonomy = read_taxonomy()
    metadata = read_metadata()
    counts = read_counts()

    # Merge counts at the indicated level
    counts, taxonomy = merge_counts(counts, taxonomy)

    # Only write out samples which are shared
    shared_samples = list(
        set(counts.columns.values) &
        set(metadata.index.values)
    )
    shared_samples.sort()
    log(f"Samples with metadata: {metadata.shape[0]:,}")
    log(f"Samples with counts: {counts.shape[1]:,}")
    log(f"Samples shared: {len(shared_samples):,}")

    assert len(shared_samples) > 1

    # Write out files
    metadata.reindex(index=shared_samples).to_csv("samplesheet.csv")
    counts = counts.reindex(columns=shared_samples).T
    counts.to_csv("counts.csv")
    counts.apply(lambda r: r / r.sum(), axis=1).to_csv("proportions.csv")
    taxonomy.to_csv("taxonomy.csv")
