#!/usr/bin/env python3

import os
import pandas as pd
import sys
from pathlib import Path
import logging


class MergeMetaphlan:

    logger: logging.Logger

    def __init__(self):
        self.setup_logger()
        self.get_args()
        self.merge_abunds()
        self.write()

    def setup_logger(self):
        self.logger = logging.getLogger()
        self.logger.setLevel(logging.INFO)
        console_handler = logging.StreamHandler()
        log_formatter = logging.Formatter(
            '%(asctime)s %(levelname)-8s [MergeMetaphlan] %(message)s'
        )
        console_handler.setFormatter(log_formatter)
        self.logger.addHandler(console_handler)

    def get_args(self):

        self.path = sys.argv[1]
        self.logger.info(f"Reading files from {self.path}")

        self.file_suffix = os.getenv("FILE_SUFFIX", ".metaphlan")
        self.logger.info(f"File suffix: {self.file_suffix}")

        self.tax_level = os.getenv("TAX_LEVEL", "species")
        self.logger.info(f"Taxonomic level: {self.tax_level}")

        self.metric = os.getenv(
            "METRIC",
            "estimated_number_of_reads_from_the_clade"
        )
        self.logger.info(f"Metric: {self.metric}")

        self.output_fp = os.getenv("OUTPUT_PATH", "merged.csv")
        self.logger.info(f"Output path: {self.output_fp}")

    def merge_abunds(self):
        self.abund = pd.DataFrame({
            name: abund
            for name, abund in self.yield_abunds()
        }).fillna(0)

        msg = f"Error: No data found in {self.path}"
        assert self.abund.shape[1] > 0, msg

    def yield_abunds(self):
        for path in Path(self.path).rglob(f"*{self.file_suffix}"):
            sample_name = path.name[:-(len(self.file_suffix))]
            self.logger.info(f"Sample: {sample_name} - File: {path}")
            yield sample_name, self.parse_abund(path)

    def parse_abund(self, path):
        # Get the headers
        header = self.get_header(path)
        self.logger.info(f"Header: {','.join(header)}")

        # Read the table of abundances
        abund = pd.read_csv(
            path,
            comment="#",
            names=header,
            sep="\t"
        )
        self.logger.info(f"Read in {abund.shape[0]:,} lines")

        # Annotate the taxonomic level
        abund = abund.assign(
            tax_level=abund["clade_name"].apply(self.get_tax_level)
        )

        # Strip out the taxonomic path from the clade name
        abund = abund.assign(
            clade_name=abund["clade_name"].apply(
                lambda s: s.rsplit("|", 1)[-1]
            )
        )

        # Filter to the specified level
        abund = abund.query(f"tax_level == '{self.tax_level}'")
        self.logger.info(f"Filtered down to {abund.shape[0]:,} lines")

        # Only return the specified metric
        msg = f"Did not find {self.metric} in header"
        assert self.metric in abund.columns.values, msg

        return abund.set_index("clade_name")[self.metric]

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
                    return line[1:].rstrip("\n").split("\t")

    def write(self):
        self.logger.info(f"Writing out to {self.output_fp}")
        self.abund.to_csv(self.output_fp)


if __name__ == "__main__":
    MergeMetaphlan()
