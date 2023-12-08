#!/usr/bin/env python3

import pandas as pd
import logging


class ReadMetaphlan:

    logger: logging.Logger
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
        self.setup_logger()
        self.get_args()
        self.read_abund()
        self.write_abund()
        self.write_taxonomy()

    def setup_logger(self):
        self.logger = logging.getLogger()
        self.logger.setLevel(logging.INFO)
        console_handler = logging.StreamHandler()
        log_formatter = logging.Formatter(
            '%(asctime)s %(levelname)-8s [ReadMetaphlan] %(message)s'
        )
        console_handler.setFormatter(log_formatter)
        self.logger.addHandler(console_handler)

    def get_args(self):

        self.path = "${file}"
        self.logger.info(f"Reading from {self.path}")

        self.name = "${sample}"
        self.logger.info(f"Name: {self.name}")

        self.tax_level = "${params.tax_level}"
        self.logger.info(f"Taxonomic level: {self.tax_level}")

        self.metric = "${params.metric}"
        self.logger.info(f"Metric: {self.metric}")

        self.output_fp = f"{self.name}.csv"
        self.logger.info(f"Output path: {self.output_fp}")

    def write_taxonomy(self):
        # Format the taxonomy table
        self.logger.info("Formatting taxonomy")
        tax = pd.DataFrame([
            self.parse_taxonomy(tax_string)
            for tax_string in self.abund["clade_name"].values
        ])
        tax = (
            tax
            .assign(index=tax['species'])
            .set_index('index')
        )
        # Write out the taxonomy table
        self.logger.info("Writing out to taxonomy.csv")
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
        self.logger.info(f"Header: {','.join(header)}")

        # Read the table of abundances
        self.abund = pd.read_csv(
            self.path,
            comment="#",
            names=header,
            sep="\\t",
            header=None
        )
        self.logger.info(f"Read in {self.abund.shape[0]:,} lines")
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
        self.logger.info(f"Filtered down to {self.abund.shape[0]:,} lines")

        # Only return the specified metric
        msg = f"Did not find {self.metric} in header"
        assert self.metric in self.abund.columns.values, msg

        self.abund = self.abund.reindex(columns=["clade_name", self.metric])

    @staticmethod
    def get_tax_level(clade_name):
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
        self.logger.info(f"Writing out to {self.output_fp}")
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
