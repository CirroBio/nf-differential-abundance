process read {
    container "${params.container}"

    input:
    file abund
    file "metadata.csv"

    output:
    path "counts.csv", emit: counts
    path "proportions.csv", emit: proportions
    path "taxonomy.csv", emit: taxonomy
    path "samplesheet.csv", emit: samplesheet

    script:
    template "read_table.py"
}


