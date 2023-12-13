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


df = pd.read_csv("corncob_results.csv")

for parameter, dat in df.groupby("parameter"):
    if parameter.endswith("(Intercept)"):
        continue
    log(f"Writing out results for {parameter}")
    log(f"Input Rows: {dat.shape[0]:,}")
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
