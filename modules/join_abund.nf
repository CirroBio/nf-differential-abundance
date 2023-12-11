process join_abund {
    container "${params.container}"
    publishDir "${params.data_output}/${params.tool}/",
        mode: 'copy',
        overwrite: true,
        saveAs: { filename -> "${params.tax_level}.${params.metric}.csv" }

    input:
    path "inputs/"

    output:
    path "${params.metric}.csv"

    script:
    template "join_abund.py"

}