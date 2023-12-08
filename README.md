# nf-differential-abundance
Nextflow workflow for differential abundance analysis

## Inputs

The format of input data can vary widely. Because the parameters
used to specify the inputs may vary by input type, this repository
will use different `main-script` entrypoints for the different
input types.

### Shared Parameters

Parameters which are used in common across all input types.

- `metadata`: DSV with sample metadata
- `sample_orientation`: Each sample's information is found in a single `row` or `column` (default: `row`)

### Tabular Abundance Data

Script: `main.nf`

Parameters:

- 

### MetaPhlAn

Script: `metaphlan.nf`

Parameters:

- `samplesheet`: URI to CSV with columns `sample` and `file`
- `tax_level`: `'species'` by default

## Glossary

- DSV: Data Separated Values (e.g. CSV or TSV)
