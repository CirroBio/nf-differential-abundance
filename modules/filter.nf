process filter {
    container "${params.container}"
    publishDir "${params.data_output}/logs/", mode: 'copy', overwrite: true, pattern: "*.log"

    input:
    path "counts.csv"
    path "taxonomy.csv"
    path "metadata.csv"

    output:
    tuple path("filtered_counts.csv"), path("filtered_taxonomy.csv"), path("filtered_metadata.csv"), emit: csv
    path "*.log", emit: log, optional: true

    script:
    template "filter.py"
}
