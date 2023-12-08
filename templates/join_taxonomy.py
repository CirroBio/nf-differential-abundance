#!/usr/bin/env python3
import pandas as pd
from pathlib import Path


def main():
    (
        pd.concat(
            [
                pd.read_csv(fp, index_col=0)
                for fp in Path("inputs").rglob("taxonomy.*.csv")
            ]
        ).drop_duplicates()
        .sort_index(axis=0)
        .to_csv("taxonomy.csv", index_label="index")
    )


if __name__ == "__main__":
    main()
