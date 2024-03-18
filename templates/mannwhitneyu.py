#!/usr/bin/env python3

import numpy as np
import pandas as pd
from scipy import stats
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


def filter_metadata(metadata):
    filter = "${params.filter}"
    if len(filter) > 0:
        logging.info(f"Filtering samples by {filter}")
        metadata = metadata.query(filter)
        logging.info(f"Filtered to {metadata.shape[0]:,} samples")
    return metadata


def filter_counts(counts, metadata):

    # Subset the counts to only include the samples in the (filtered) metadata
    counts = counts.loc[metadata.index]

    min_prevalence = float("${params.min_prevalence}")
    logging.info(f"Minimum prevalence filter: {min_prevalence} (proportion of samples)")
    counts = counts.loc[:, (counts > 0).mean() > min_prevalence]
    logging.info(f"Filtered to {counts.shape[1]:,} taxa")
    return counts


def clr(r: pd.Series):
    logvals = np.log10(r + 1)
    return logvals - np.mean(logvals)


def mannwhitneyu(
    org: str,
    org_props: pd.Series,
    metadata: pd.DataFrame,
    formula_elem: str
):
    msg = f"{formula_elem} not found in metadata columns"
    assert formula_elem in metadata.columns.values, msg

    metadata_values = metadata[formula_elem].dropna()
    vc = ", ".join([
        f"{kw} (n={n:,})"
        for kw, n in metadata_values.value_counts().items()
    ])
    msg = f"mannwhitneyu test only appropriate for two categories, not {vc}"
    assert metadata_values.unique().shape[0] == 2, msg

    grouped_values = [
        vals
        for _, vals in org_props.groupby(metadata_values)
    ]

    res = stats.mannwhitneyu(grouped_values[0], grouped_values[1])

    # Calculate the log-fold-change
    lfc = np.mean(grouped_values[1]) - np.mean(grouped_values[0])

    return dict(
        org=org,
        p_value=res.pvalue,
        statistic=res.statistic,
        lfc=lfc,
        metadata=formula_elem
    )


if __name__ == "__main__":

    metadata = pd.read_csv("metadata.csv", index_col=0)
    log(metadata.head().to_csv())

    log("Filtering metadata")
    metadata = filter_metadata(metadata)

    counts = pd.read_csv("counts.csv", index_col=0)
    log(counts.head().to_csv())

    counts = counts.loc[metadata.index]

    log("Filtering counts")
    counts = filter_counts(counts, metadata)

    log("Calculating CLR")

    props = counts.apply(clr, axis=1)
    log(props.head().to_csv())

    results = pd.DataFrame([
        mannwhitneyu(
            org,
            org_props,
            metadata,
            formula_elem.strip()
        )
        for formula_elem in "${params.formula}".split("+")
        for org, org_props in props.items()
    ])
    results.set_index("org", inplace=True)

    log(results.head().to_csv())

    for metadata, df in results.groupby("metadata"):
        df.to_csv(f"{metadata}.csv")
