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


df = pd.read_csv("radEmu_results.csv")

for parameter, dat in df.groupby("covariate"):
    log(f"Writing out results for {parameter}")
    log(f"Input Rows: {dat.shape[0]:,}")
    (
        dat
        .rename(columns={
            "pval": "p_value",
            "se": "std_error",
            "estimate": "est_coef"
        })
        .drop(
            columns=["covariate", "category_num"]
        )
        .sort_values(by="p_value")
        .to_csv(f"{parameter}.csv", index=None)
    )
