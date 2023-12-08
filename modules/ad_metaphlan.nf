process ad_metaphlan {
    container "${params.container}"
    publishDir "${params.output}/anndata/", mode: 'copy', overwrite: true, pattern: "*.h5ad"

    input:
        path "counts.csv"
        path "proportions.csv"
        path "taxonomy.csv"
        path "samplesheet.csv"
        path "corncob/"

    output:
        tuple path("metaphlan.h5ad", "variables.txt")

    script:
        template "ad_metaphlan.py"
}