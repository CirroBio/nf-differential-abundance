#!/usr/bin/env python3

import json
from typing import Tuple
import pandas as pd
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


def read_counts() -> Tuple[pd.DataFrame, pd.DataFrame]:

    log("Parsing read alignments")
    df = pd.read_csv("read_alignments.csv.gz")
    log(df.head().to_csv())
    log(f"Number of lines: {df.shape[0]:,}")

    # Filter by minimum coverage
    log("Filtering by minimum coverage: ${params.min_coverage}")
    df = df.query("coverage >= ${params.min_coverage}")
    log(f"Number of lines: {df.shape[0]:,}")

    # Get a wide table of read counts
    wide_df = (
        df
        .pivot_table(
            index="specimen",
            columns="id",
            values="nreads"
        )
        .fillna(0)
    )
    log(f"Number of genes: {wide_df.shape[1]:,}")
    log(f"Number of samples: {wide_df.shape[0]:,}")

    # If unaligned reads should be included in the analysis
    if "${params.use_unaligned}".lower() == "true":
        log("Including unaligned reads in the analysis")

        # Get the total number of reads per specimen
        tot_reads = (
            df
            .reindex(columns=["specimen", "tot_reads"])
            .drop_duplicates()
            .set_index("specimen")
            ["tot_reads"]
        )
        log(tot_reads.to_csv())

        # Add a column for the number of unaligned reads
        wide_df = wide_df.assign(
            UNALIGNED=tot_reads-wide_df.sum(axis=1)
        )
        log(wide_df["UNALIGNED"].to_csv())

    # If only aligned reads should be used
    else:
        log("Only including aligned reads in the analysis")

    # Calculate the proportion of reads which align to each
    log("Calculating proportional abundance of all genes")
    prop_df = (wide_df.T / wide_df.sum(axis=1)).T

    return wide_df, prop_df


def read_metadata():

    # Reading the sample metadata
    log("Reading the provided metadata")
    metadata = pd.read_csv("metadata.csv")
    log(metadata.head().to_csv())
    msg = "Could not find 'sample' column"
    assert "sample" in metadata.columns.values, msg
    log("Setting index to the 'sample' column")
    metadata = metadata.set_index("sample")

    return metadata


def main():

    metadata = read_metadata()
    counts, proportions = read_counts()

    # Only write out samples which are shared
    shared_samples = list(
        set(counts.index.values) &
        set(metadata.index.values)
    )
    shared_samples.sort()
    log(f"Samples with metadata: {metadata.shape[0]:,}")
    log(f"Samples with counts: {counts.shape[0]:,}")
    log(f"Samples shared: {len(shared_samples):,}")

    assert len(shared_samples) > 1

    metadata = metadata.reindex(index=shared_samples)
    counts = counts.reindex(index=shared_samples)
    proportions = proportions.reindex(index=shared_samples)

    # Format a dummy taxonomy
    taxonomy = format_taxonomy(counts.columns.values)

    log("Writing out the sample metadata")
    metadata.to_csv("samplesheet.csv")

    for df, kw in [
        (counts, "counts"),
        (proportions, "proportions"),
        (taxonomy, "taxonomy")
    ]:
        log(f"Writing out {kw} table")
        log(df.head().to_csv())
        df.to_csv(f"{kw}.csv")


def format_taxonomy(gene_ids: list):

    tax_levels = [
        "kingdom",
        "phylum",
        "class",
        "order",
        "family",
        "genus",
        "species",
        "strain"
    ]

    return (
        pd.DataFrame([
            dict(
                name=gene_id,
                **{
                    level: gene_id
                    for level in tax_levels
                }
            )
            for gene_id in gene_ids
        ])
        .set_index("name")
        .reindex(columns=tax_levels)
    )


if __name__ == "__main__":
    main()
