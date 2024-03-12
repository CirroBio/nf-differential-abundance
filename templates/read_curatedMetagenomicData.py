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


def split_data(df: pd.DataFrame) -> Tuple[pd.DataFrame, pd.DataFrame]:
    abund_cnames = [
        cname
        for cname in df.columns.values
        if "__" in cname
    ]
    metadata_cnames = [
        cname for cname in df.columns.values
        if cname not in abund_cnames
    ]
    metadata = df.reindex(columns=metadata_cnames)
    abunds = df.reindex(columns=abund_cnames)

    # Only keep samples with abunds > 0
    abunds = abunds.loc[abunds.sum(axis=1) > 0]
    metadata = metadata.loc[abunds.index]

    return metadata, abunds


def is_counts(abunds) -> bool:
    return (
        abunds.applymap(
            lambda v: isinstance(v, int)
        )
        .all()
        .all()
    )


def main():

    log("Reading in full table")
    df = pd.read_csv("input.csv", sep="\t", index_col=0)
    log(df.head().to_csv())

    log("Splitting metadata and abundances")
    metadata, abunds = split_data(df)

    log("Writing out the sample metadata")
    metadata.to_csv("samplesheet.csv")

    abunds, taxonomy = filter_tax_level(abunds)

    log("Making counts and proportions")
    counts, proportions = make_counts_proportions(abunds)

    for df, kw in [
        (counts, "counts"),
        (proportions, "proportions"),
        (taxonomy, "taxonomy")
    ]:
        log(f"Writing out {kw} table")
        log(df.head().to_csv())
        df.to_csv(f"{kw}.csv")


def filter_tax_level(abunds: pd.DataFrame) -> Tuple[pd.DataFrame, pd.DataFrame]:

    tax_level = "${params.tax_level}"

    # Make a table with the taxonomy for each organism
    log("Parsing taxonomy")
    taxonomy = pd.DataFrame([
        parse_tax_string(tax_string)
        for tax_string in abunds.columns.values
    ]).set_index("lineage")
    log(taxonomy.head().to_csv())
    log(json.dumps(taxonomy["level"].value_counts().to_dict(), indent=4))

    # Filter to the indicated level
    log(f"Filtering to {tax_level} level")
    taxonomy = taxonomy.query(f"level == '{tax_level}'")
    log(taxonomy.head().to_csv())
    assert taxonomy.shape[0] > 0, "No organisms detected at this level"

    # Filter the abundance information as well
    abunds = abunds.reindex(columns=taxonomy.index.values)
    log("Filtered abundances")
    log(abunds.head().to_csv())

    # Rename the abundances for the organism name
    abunds = abunds.rename(
        columns=taxonomy[tax_level]
    )
    log("Renamed abundances")
    log(abunds.head().to_csv())

    # Rename the taxonomy for the specified level
    taxonomy = (
        taxonomy
        .reset_index(drop=True)
        .set_index("name")
    )
    log("Reformatted the taxonomy")
    log(taxonomy.head().to_csv())

    return abunds, taxonomy


def parse_tax_string(lineage: str) -> dict:

    mapping = dict(
        k="kingdom",
        p="phylum",
        c="class",
        o="order",
        f="family",
        g="genus",
        s="species",
        t="strain"
    )

    # Loop over each of the specified levels
    dat = dict(
        lineage=lineage
    )

    for level in lineage.split("|"):

        # Get the first character
        # slc = single-letter code
        slc = level[0]

        assert slc in mapping, f"Couldn't parse tax_level: {lineage}"

        msg = f"Expected '__' delim in '{level}': {lineage}"
        assert level[1:3] == "__", msg

        dat[mapping[slc]] = level[3:]

    # Add the information showing the final level
    dat["level"] = mapping[slc]
    dat["name"] = level[3:]
    return dat


def make_counts_proportions(abunds) -> Tuple[pd.DataFrame, pd.DataFrame]:
    if is_counts(abunds):
        log("Input data appears to be counts")
        log("Calculating proportions")
        proportions = abunds.apply(lambda r: r / r.sum(), axis=1)
        counts = abunds
        log(proportions.head().to_csv())
    else:
        log("Input data appears to be proportions")
        log("Estimating counts")
        counts = abunds.apply(
            lambda r: (r / r[r > 0].min()).apply(int), axis=1
        )

        log(counts.head().to_csv())
        proportions = abunds.apply(
            lambda r: r / r.sum(), axis=1
        )

    return counts, proportions


if __name__ == "__main__":
    main()
