process make_anndata {
    container "${params.container}"
    publishDir "${params.data_output}/anndata/", mode: 'copy', overwrite: true, pattern: "*.h5ad"

    input:
        path "counts.csv"
        path "proportions.csv"
        path "taxonomy.csv"
        path "samplesheet.csv"
        path "stats/"

    output:
        tuple path("${params.tool}.h5ad"), path("output_config.json")

    script:
        template "make_anndata.py"
}