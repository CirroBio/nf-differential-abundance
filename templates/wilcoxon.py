#!/usr/bin/env python3

import numpy as np
import pandas as pd
from scipy import stats

counts = pd.read_csv("counts.csv", index_col=0)
print(counts)

print("Calculating CLR")


def clr(r: pd.Series):
    logvals = np.log10(r + 1)
    return logvals - np.mean(logvals)


props = counts.apply(clr, axis=1)
print(props)

metadata = pd.read_csv("metadata.csv", index_col=0)
print(metadata)


def wilcoxon(
    org: str,
    org_props: pd.Series,
    metadata: pd.DataFrame,
    formula_elem: str
):
    msg = f"{formula_elem} not found in metadata columns"
    assert formula_elem in metadata.columns.values, msg

    metadata_values = metadata[formula_elem].dropna()
    vc = ", ".join([
        f"{kw} (n={n:,})"
        for kw, n in metadata_values.value_counts().items()
    ])
    msg = f"Wilcoxon test only appropriate for two categories, not {vc}"
    assert metadata_values.unique().shape[0] == 2, msg

    grouped_values = [
        vals
        for _, vals in org_props.groupby(metadata_values)
    ]

    res = stats.wilcoxon(grouped_values[0], grouped_values[1])

    # Calculate the log-fold-change
    lfc = np.mean(grouped_values[1]) - np.mean(grouped_values[0])

    return dict(
        org=org,
        p_value=res.pvalue,
        statistic=res.statistic,
        lfc=lfc,
        metadata=formula_elem
    )


results = pd.DataFrame([
    wilcoxon(
        org,
        org_props,
        metadata,
        formula_elem.strip()
    )
    for formula_elem in "${params.formula}".split("+")
    for org, org_props in props.items()
])
results.set_index("org", inplace=True)

print(results)

for metadata, df in results.groupby("metadata"):
    df.to_csv(f"{metadata}.csv")
