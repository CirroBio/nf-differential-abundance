process join_abund {
    container "${params.container}"
    publishDir "${params.data_output}/${params.tool}/",
        pattern: "*.csv",
        mode: 'copy',
        overwrite: true,
        saveAs: { filename -> "${params.tax_level}.${params.metric}.csv" }
    publishDir "${params.data_output}/logs/", mode: 'copy', overwrite: true, pattern: "*.log"

    input:
    path "inputs/"

    output:
    path "${params.metric}.csv", emit: csv
    path "*.log", emit: log, optional: true

    script:
    template "join_abund.py"

}