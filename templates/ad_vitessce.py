#!/usr/bin/env python3

import json
import anndata as ad

adata = ad.read_h5ad("adata.h5ad")
# Write out in both orientations
adata.write_zarr("obs.zarr")
adata.T.write_zarr("var.zarr")


# Format each of the visualizations
def make_vt_config():

    return {
        "version": "1.0.16",
        "name": "Microbiome ~ CRC",
        "description": "Summary of organisms associated with CRC",
        "datasets": [
            {
                "uid": "A",
                "name": "Species Level Abundances (observations)",
                "files": [
                    {
                        "fileType": "anndata.zarr",
                        "url": "http://localhost:9000/obs.zarr",
                        "options": {
                            "obsEmbedding": [
                                {
                                    "path": "obsm/X_umap",
                                    "dims": [
                                        0,
                                        1
                                    ],
                                    "embeddingType": "UMAP"
                                }
                            ],
                            "obsSets": [
                                {
                                    "name": "CRC",
                                    "path": "obs/CRC"
                                },
                                {
                                    "name": "BioProject",
                                    "path": "obs/BioProject"
                                }
                            ],
                            "obsFeatureMatrix": {
                                "path": "X"
                            }
                        },
                        "coordinationValues": {
                            "obsType": "sample",
                            "featureType": "species"
                        }
                    }
                ]
            },
            {
                "uid": "B",
                "name": "Species Level Abundances (features)",
                "files": [
                    {
                        "fileType": "anndata.zarr",
                        "url": "http://localhost:9000/var.zarr",
                        "options": {
                            "obsEmbedding": [
                                {
                                    "path": "obsm/CRC_volcano",
                                    "dims": [
                                        0,
                                        1
                                    ],
                                    "embeddingType": "CRC - Volcano Plot"
                                },
                                {
                                    "path": "obsm/CRC_ma",
                                    "dims": [
                                        0,
                                        1
                                    ],
                                    "embeddingType": "CRC - MA Plot"
                                }
                            ],
                            "obsSets": [
                                {
                                    "name": "CRC-Associated",
                                    "path": "obsm/mu.CRC/sig_diff"
                                },
                                {
                                    "name": "Mean Proportion (%)",
                                    "path": "obs/mean_proportion"
                                }
                            ],
                            "obsFeatureMatrix": {
                                "path": "X"
                            }
                        },
                        "coordinationValues": {
                            "obsType": "species",
                            "featureType": "sample"
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
                "A": "UMAP",
                "B": "CRC - Volcano Plot",
                "C": "CRC - MA Plot"
            },
            "embeddingZoom": {
                "A": 3,
                "B": 5,
                "C": 5
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
                "A": ["Escherichia_coli"],
                "B": ["None"],
                "C": ["None"]
            },
            "obsSetColor": {
                "A": null,
                "B": null
            },
            "obsSetSelection": {
                "A": null,
                "B": null
            },
            "obsType": {
                "A": "sample",
                "B": "species"
            },
            "featureType": {
                "A": "species",
                "B": "sample"
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
with open("config.json", "r") as handle:
    config: dict = json.load(handle)

for kw, elem in config.items():
    with open(f"{kw}.vt.json", "w") as handle:
        json.dump(elem, handle, indent=4)
