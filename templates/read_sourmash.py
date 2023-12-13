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


def log(s: str):
    for line in s.split("\\n"):
        logging.info(line)


class ReadSourmash:

    levels = [
        'kingdom',
        'phylum',
        'class',
        'order',
        'family',
        'genus',
        'species'
    ]

    def __init__(self):
        self.get_args()
        self.read_abund()
        self.write_abund("bp_match_at_rank", "counts")
        self.write_abund("f_weighted_at_rank", "proportions")
        self.write_taxonomy()

    def get_args(self):

        self.path = "${file}"
        log(f"Reading from {self.path}")

        self.name = "${sample}"
        log(f"Name: {self.name}")

        self.tax_level = "${params.tax_level}"
        log(f"Taxonomic level: {self.tax_level}")

    def write_taxonomy(self):
        # Format the taxonomy table
        log("Formatting taxonomy")
        tax = pd.DataFrame([
            self.parse_taxonomy(tax_string)
            for tax_string in self.abund["lineage"].values
        ])
        tax = (
            tax
            .assign(index=tax[self.tax_level])
            .set_index('index')
            .rename(
                index=lambda s: s.replace(" ", "_")
            )
        )
        # Write out the taxonomy table
        log("Writing out to taxonomy.csv")
        tax.to_csv("taxonomy.csv")

    def parse_taxonomy(self, tax_string):
        if not isinstance(tax_string, str):
            print(tax_string)
        anc = {
            self.get_tax_level(taxon): taxon[3:]
            for taxon in tax_string.split(";")
        }
        # Fill in any missing values
        for i, n in enumerate(self.levels):
            if anc.get(n) is None:
                anc[n] = self.levels[i - 1]
        return anc

    def read_abund(self):

        # Read the table of abundances
        self.abund = pd.read_csv(self.path)
        log(f"Read in {self.abund.shape[0]:,} lines")

        # Filter to the specified level
        self.abund = self.abund.query(f"rank == '{self.tax_level}'")
        log(f"Filtered down to {self.abund.shape[0]:,} lines")

        # Just get the columns of interest
        self.abund = self.abund.reindex(
            columns=["lineage", "f_weighted_at_rank", "bp_match_at_rank"]
        )

        # Remove any unclassified
        self.abund = self.abund.query("lineage != 'unclassified'")

    @staticmethod
    def get_tax_level(lineage: str):
        if lineage.upper() == "UNCLASSIFIED":
            return

        mapping = dict(
            d="kingdom",
            p="phylum",
            c="class",
            o="order",
            f="family",
            g="genus",
            s="species",
            t="strain"
        )
        # Get the first character of the final label
        # slc = single-letter code
        slc = lineage[0]

        assert slc in mapping, f"Couldn't parse tax_level: {lineage}"

        return mapping[slc]

    def write_abund(self, metric, suffix):
        output_fp = f"{self.name}.{suffix}.csv"
        log(f"Writing out to {output_fp}")
        # Only write out the final organism name
        (
            self.abund
            .assign(
                lineage=self.abund["lineage"].apply(
                    lambda s: s.split(";")[-1][3:].replace(" ", "_")
                )
            )
            .rename(
                columns={
                    metric: self.name
                }
            )
            .reindex(columns=["lineage", self.name])
            .to_csv(
                output_fp,
                index=None
            )
        )


if __name__ == "__main__":
    ReadSourmash()
