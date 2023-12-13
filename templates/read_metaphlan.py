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


class ReadMetaphlan:

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
        self.write_abund()
        self.write_taxonomy()

    def get_args(self):

        self.path = "${file}"
        log(f"Reading from {self.path}")

        self.name = "${sample}"
        log(f"Name: {self.name}")

        self.tax_level = "${params.tax_level}"
        log(f"Taxonomic level: {self.tax_level}")

        self.metric = "${params.metric}"
        log(f"Metric: {self.metric}")

        self.output_fp = f"{self.name}.csv"
        log(f"Output path: {self.output_fp}")

    def write_taxonomy(self):
        # Format the taxonomy table
        log("Formatting taxonomy")
        tax = pd.DataFrame([
            self.parse_taxonomy(tax_string)
            for tax_string in self.abund["clade_name"].values
        ])
        tax = (
            tax
            .assign(index=tax[self.tax_level])
            .set_index('index')
        )
        # Write out the taxonomy table
        log("Writing out to taxonomy.csv")
        tax.to_csv("taxonomy.csv")

    def parse_taxonomy(self, tax_string):
        if not isinstance(tax_string, str):
            print(tax_string)
        anc = {
            self.get_tax_level(taxon[0]): taxon[3:]
            for taxon in tax_string.split("|")
        }
        # Fill in any missing values
        for i, n in enumerate(self.levels):
            if anc.get(n) is None:
                anc[n] = self.levels[i - 1]
        return anc

    def read_abund(self):
        # Get the headers
        header = self.get_header(self.path)
        log(f"Header: {','.join(header)}")

        # Read the table of abundances
        self.abund = pd.read_csv(
            self.path,
            comment="#",
            names=header,
            sep="\\t",
            header=None
        )
        log(f"Read in {self.abund.shape[0]:,} lines")
        assert (
            self.abund["clade_name"]
            .apply(lambda v: isinstance(v, str))
            .all()
        ), self.abund["clade_name"]

        # Annotate the taxonomic level
        self.abund = self.abund.assign(
            tax_level=self.abund["clade_name"].apply(self.get_tax_level)
        )

        # Filter to the specified level
        self.abund = self.abund.query(f"tax_level == '{self.tax_level}'")
        log(f"Filtered down to {self.abund.shape[0]:,} lines")

        # Only return the specified metric
        msg = f"Did not find {self.metric} in header"
        assert self.metric in self.abund.columns.values, msg

        self.abund = self.abund.reindex(columns=["clade_name", self.metric])

    @staticmethod
    def get_tax_level(clade_name):
        if clade_name == "UNCLASSIFIED":
            return

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
        # Get the first character of the final label
        # slc = single-letter code
        slc = clade_name.split("|")[-1][0]

        assert slc in mapping, f"Couldn't parse tax_level: {clade_name}"

        return mapping[slc]

    @staticmethod
    def get_header(path):
        with open(path) as handle:
            for line in handle:
                if line.startswith("#clade_name"):
                    return line[1:].rstrip("\\n").split("\\t")

    def write_abund(self):
        # Only write out the final organism name
        log(f"Writing out to {self.output_fp}")
        (
            self.abund
            .assign(
                clade_name=self.abund["clade_name"].apply(
                    lambda s: s.split("|")[-1][3:]
                )
            )
            .rename(
                columns={
                    self.metric: self.name
                }
            )
            .to_csv(
                self.output_fp,
                index=None
            )
        )


if __name__ == "__main__":
    ReadMetaphlan()
