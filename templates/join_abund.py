#!/usr/bin/env python3
import pandas as pd
from pathlib import Path


def read(fp):

    df = pd.read_csv(fp, index_col=0)
    for cname, cvals in df.items():
        yield cname, cvals


def main():
    (
        pd.DataFrame(
            {
                name: vals
                for fp in Path("inputs").rglob("*")
                for name, vals in read(fp)
            }
        )
        .T
        .fillna(0)
        .sort_index(axis=0)
        .sort_index(axis=1)
        .to_csv("abund.csv")
    )


if __name__ == "__main__":
    main()
