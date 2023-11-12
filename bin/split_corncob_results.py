#!/usr/bin/env python3

import pandas as pd
import sys


def main():
    input_fp = sys.argv[1]
    all_df = pd.read_csv(input_fp)
    for parameter, df in all_df.groupby("parameter"):
        output_fp = f"{parameter.replace('(').replace(')')}.csv"
        (
            df
            .pivot(
                index="id",
                columns="variable",
                values="value"
            )
            .dropna()
            .to_csv(output_fp)
        )


if __name__ == "__main__":
    main()
