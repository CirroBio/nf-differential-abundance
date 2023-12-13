#!/usr/bin/env python3

import json
import anndata as ad

adata = ad.read_h5ad("adata.h5ad")
# Write out in both orientations
adata.write_zarr("obs.zarr")
adata.T.write_zarr("var.zarr")


# Format each of the visualizations
def make_vt_config(kwargs):

    return {
        "version": "1.0.16",
        "name": kwargs["title"],
        "description": kwargs["description"],
        "datasets": [
            {
                "uid": "A",
                "name": f"{kwargs['obs_title']} ({kwargs['obs_type'].title()})",
                "files": [
                    {
                        "fileType": "anndata.zarr",
                        "url": "obs.zarr",
                        "options": {
                            "obsEmbedding": [
                                {
                                    "path": kwargs["ordination_path"],
                                    "dims": [
                                        0,
                                        1
                                    ],
                                    "embeddingType": kwargs["ordination_title"]
                                }
                            ],
                            "obsSets": kwargs["obs_sets"],
                            "obsFeatureMatrix": {
                                "path": "X"
                            }
                        },
                        "coordinationValues": {
                            "obsType": kwargs["obs_type"],
                            "featureType": kwargs["feature_type"]
                        }
                    }
                ]
            },
            {
                "uid": "B",
                "name": f"{kwargs['obs_title']} ({kwargs['feature_type'].title()})",
                "files": [
                    {
                        "fileType": "anndata.zarr",
                        "url": "var.zarr",
                        "options": {
                            "obsEmbedding": [
                                {
                                    "path": f"obsm/{kwargs['volcano_varm']}",
                                    "dims": [
                                        0,
                                        1
                                    ],
                                    "embeddingType": kwargs["volcano_title"]
                                },
                                {
                                    "path": f"obsm/{kwargs['ma_varm']}",
                                    "dims": [
                                        0,
                                        1
                                    ],
                                    "embeddingType": kwargs["ma_title"]
                                }
                            ],
                            "obsSets": kwargs["var_sets"],
                            "obsFeatureMatrix": {
                                "path": "X"
                            }
                        },
                        "coordinationValues": {
                            "obsType": kwargs["feature_type"],
                            "featureType": kwargs["obs_type"]
                        }
                    }
                ]
            }
        ],
        "coordinationSpace": {
            "dataset": {
                "A": "A",
                "B": "B"
            },
            "embeddingType": {
                "A": kwargs["ordination_title"],
                "B": kwargs["volcano_title"],
                "C": kwargs["ma_title"]
            },
            "embeddingZoom": {
                "A": 8,
                "B": 8,
                "C": 8
            },
            "embeddingTargetX": {
                "A": 0,
                "B": 0,
                "C": 0
            },
            "embeddingTargetY": {
                "A": 0,
                "B": 0,
                "C": 0
            },
            "embeddingObsRadiusMode": {
                "A": "manual",
                "B": "manual",
                "C": "manual"
            },
            "embeddingObsRadius": {
                "A": 2,
                "B": 2,
                "C": 2
            },
            "obsColorEncoding": {
                "A": "cellSetSelection",
                "B": "geneSelection",
                "C": "cellSetSelection"
            },
            "featureSelection": {
                "A": [kwargs["top_taxon"]],
                "B": ["None"],
                "C": ["None"]
            },
            "obsSetColor": {
                "A": None,
                "B": None
            },
            "obsSetSelection": {
                "A": None,
                "B": None
            },
            "obsType": {
                "A": kwargs["obs_type"],
                "B": kwargs["feature_type"]
            },
            "featureType": {
                "A": kwargs["feature_type"],
                "B": kwargs["obs_type"]
            }
        },
        "layout": [
            {
                "component": "heatmap",
                "coordinationScopes": {
                    "dataset": "A",
                    "obsType": "A",
                    "featureType": "A",
                    "featureSelection": "A",
                    "obsSetColor": "A",
                    "obsSetSelection": "A"
                },
                "x": 0,
                "y": 0,
                "w": 4,
                "h": 6
            },
            {
                "component": "scatterplot",
                "coordinationScopes": {
                    "dataset": "A",
                    "obsType": "A",
                    "featureType": "A",
                    "embeddingZoom": "A",
                    "embeddingTargetX": "A",
                    "embeddingTargetY": "A",
                    "embeddingType": "A",
                    "embeddingObsRadiusMode": "A",
                    "embeddingObsRadius": "A",
                    "obsColorEncoding": "A",
                    "featureSelection": "B",
                    "obsSetColor": "A",
                    "obsSetSelection": "A"
                },
                "x": 4,
                "y": 0,
                "w": 3,
                "h": 6
            },
            {
                "component": "scatterplot",
                "coordinationScopes": {
                    "dataset": "A",
                    "obsType": "A",
                    "featureType": "A",
                    "embeddingZoom": "A",
                    "embeddingTargetX": "A",
                    "embeddingTargetY": "A",
                    "embeddingType": "A",
                    "embeddingObsRadiusMode": "A",
                    "embeddingObsRadius": "A",
                    "obsColorEncoding": "B",
                    "featureSelection": "A",
                    "obsSetColor": "A",
                    "obsSetSelection": "A"
                },
                "x": 7,
                "y": 0,
                "w": 3,
                "h": 6
            },
            {
                "component": "obsSetSizes",
                "coordinationScopes": {
                    "dataset": "A",
                    "obsType": "A",
                    "featureType": "A",
                    "obsSetColor": "A",
                    "obsSetSelection": "A"
                },
                "x": 10,
                "y": 0,
                "w": 2,
                "h": 6
            },
            {
                "component": "scatterplot",
                "coordinationScopes": {
                    "dataset": "B",
                    "obsType": "B",
                    "featureType": "B",
                    "embeddingZoom": "B",
                    "embeddingTargetX": "B",
                    "embeddingTargetY": "B",
                    "embeddingType": "B",
                    "embeddingObsRadiusMode": "B",
                    "embeddingObsRadius": "B",
                    "obsColorEncoding": "C",
                    "featureSelection": "C",
                    "obsSetColor": "B",
                    "obsSetSelection": "B"
                },
                "x": 0,
                "y": 6,
                "w": 4,
                "h": 6
            },
            {
                "component": "scatterplot",
                "coordinationScopes": {
                    "dataset": "B",
                    "obsType": "B",
                    "featureType": "B",
                    "embeddingZoom": "C",
                    "embeddingTargetX": "C",
                    "embeddingTargetY": "C",
                    "embeddingType": "C",
                    "embeddingObsRadiusMode": "C",
                    "embeddingObsRadius": "C",
                    "obsColorEncoding": "C",
                    "featureSelection": "C",
                    "obsSetColor": "B",
                    "obsSetSelection": "B"
                },
                "x": 4,
                "y": 6,
                "w": 3,
                "h": 6
            },
            {
                "component": "obsSetFeatureValueDistribution",
                "coordinationScopes": {
                    "dataset": "A",
                    "obsType": "A",
                    "featureType": "A",
                    "featureSelection": "A",
                    "obsSetColor": "A",
                    "obsSetSelection": "A"
                },
                "x": 7,
                "y": 6,
                "w": 3,
                "h": 6
            },
            {
                "component": "featureList",
                "coordinationScopes": {
                    "dataset": "A",
                    "obsType": "A",
                    "featureType": "A",
                    "featureSelection": "A"
                },
                "x": 10,
                "y": 6,
                "w": 2,
                "h": 6
            }
        ],
        "initStrategy": "auto"
    }


# Read in each of the configurations
with open("input_config.json", "r") as handle:
    config: dict = json.load(handle)

print(config)

# Keep a list of all of the files which have been written
chart_manifest = []
for kw, elem in config.items():
    vt_config = make_vt_config(elem)
    fp = f"{kw}.vt.json"
    with open(fp, "w") as handle:
        json.dump(vt_config, handle, indent=4)

    chart_manifest.append(
        dict(
            type="vitessce",
            config=fp,
            name=elem["title"],
            desc=elem["description"]
        )
    )

with open("chart.manifest.json", "w") as handle:
    json.dump(chart_manifest, handle, indent=4)
