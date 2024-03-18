#!/usr/bin/env python3

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


def main():
    counts = read_table("counts.csv")
    metadata = read_table("metadata.csv")
    taxonomy = read_table("taxonomy.csv")

    # Filter the metadata, if a filter was provided by the user
    metadata = filter_metadata(metadata)

    # Filter the counts, if min_prevalence > 0
    counts = filter_counts(counts, metadata)

    # Filter the taxonomy to match the filtered counts
    taxonomy = taxonomy.loc[counts.columns]

    # Write out the filtered tables
    counts.to_csv("filtered_counts.csv")
    metadata.to_csv("filtered_metadata.csv")
    taxonomy.to_csv("filtered_taxonomy.csv")


def filter_counts(counts, metadata):

    # Subset the counts to only include the samples in the (filtered) metadata
    counts = counts.loc[metadata.index]

    min_prevalence = float("${params.min_prevalence}")
    logging.info(f"Minimum prevalence filter: {min_prevalence} (proportion of samples)")
    counts = counts.loc[:, (counts > 0).mean() > min_prevalence]
    logging.info(f"Filtered to {counts.shape[1]:,} taxa")
    return counts


def filter_metadata(metadata):
    filter = "${params.filter}"
    if len(filter) > 0:
        logging.info(f"Filtering samples by {filter}")
        metadata = metadata.query(filter)
        logging.info(f"Filtered to {metadata.shape[0]:,} samples")
    return metadata


def read_table(fp):
    logging.info(f"Reading {fp}")
    df = pd.read_csv(fp, index_col=0)
    logging.info(f"Read {df.shape[0]:,} rows and {df.shape[1]:,} columns")
    logging.info(df.head().T.head().T)
    return df


if __name__ == "__main__":
    main()
