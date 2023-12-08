#!/usr/bin/env python3

import pandas as pd

df = pd.read_csv("corncob_results.csv")

for parameter, dat in df.groupby("parameter"):
    if parameter.endswith("(Intercept)"):
        continue
    (
        dat
        .pivot(
            index="id",
            columns="variable",
            values="value"
        )
        .rename(columns={
            "Pr(>|t|)": "p_value",
            "Std. Error": "std_error",
            "Estimate": "est_coef",
            "t value": "t_value"
        })
        .sort_values(by="p_value")
        .to_csv(f"{parameter}.csv")
    )
