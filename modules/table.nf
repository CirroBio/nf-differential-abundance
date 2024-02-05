process read {
    container "${params.container}"
    publishDir "${params.data_output}/logs/", mode: 'copy', overwrite: true, pattern: "*.log"

    input:
    file abund
    file "metadata.csv"

    output:
    path "counts.csv", emit: counts
    path "proportions.csv", emit: proportions
    path "taxonomy.csv", emit: taxonomy
    path "samplesheet.csv", emit: samplesheet
    path "*.log", emit: log, optional: true

    script:
    template "read_table.py"
}


