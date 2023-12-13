#!/usr/bin/env python3

from typing import Tuple
import pandas as pd


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

    print("Reading in full table")
    df = pd.read_csv("input.csv", sep="\t", index_col=0)
    print(df)

    print("Splitting metadata and abundances")
    metadata, abunds = split_data(df)

    print("Writing out the sample metadata")
    metadata.to_csv("samplesheet.csv")

    abunds, taxonomy = filter_tax_level(abunds)

    print("Making counts and proportions")
    counts, proportions = make_counts_proportions(abunds)

    for df, kw in [
        (counts, "counts"),
        (proportions, "proportions"),
        (taxonomy, "taxonomy")
    ]:
        print(f"Writing out {kw} table")
        print(df)
        df.to_csv(f"{kw}.csv")


def filter_tax_level(abunds: pd.DataFrame) -> Tuple[pd.DataFrame, pd.DataFrame]:

    tax_level = "${params.tax_level}"

    # Make a table with the taxonomy for each organism
    print("Parsing taxonomy")
    taxonomy = pd.DataFrame([
        parse_tax_string(tax_string)
        for tax_string in abunds.columns.values
    ]).set_index("lineage")
    print(taxonomy)
    print(taxonomy["level"].value_counts())

    # Filter to the indicated level
    print(f"Filtering to {tax_level} level")
    taxonomy = taxonomy.query(f"level == '{tax_level}'")
    print(taxonomy)
    assert taxonomy.shape[0] > 0, "No organisms detected at this level"

    # Filter the abundance information as well
    abunds = abunds.reindex(columns=taxonomy.index.values)
    print("Filtered abundances")
    print(abunds)

    # Rename the abundances for the organism name
    abunds = abunds.rename(
        columns=taxonomy[tax_level]
    )
    print("Renamed abundances")
    print(abunds)

    # Rename the taxonomy for the specified level
    taxonomy = (
        taxonomy
        .reset_index(drop=True)
        .set_index("name")
    )
    print("Reformatted the taxonomy")
    print(taxonomy)

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
        print("Input data appears to be counts")
        print("Calculating proportions")
        proportions = abunds.apply(lambda r: r / r.sum(), axis=1)
        counts = abunds
        print(proportions)
    else:
        print("Input data appears to be proportions")
        print("Estimating counts")
        counts = abunds.apply(
            lambda r: (r / r[r > 0].min()).apply(int), axis=1
        )

        print(counts)
        proportions = abunds.apply(
            lambda r: r / r.sum(), axis=1
        )

    return counts, proportions


if __name__ == "__main__":
    main()
