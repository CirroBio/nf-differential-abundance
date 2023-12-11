process read {
    container "${params.container}"

    input:
    file "input.csv"

    output:
    path "counts.csv", emit: counts
    path "proportions.csv", emit: proportions
    path "taxonomy.csv", emit: taxonomy
    path "samplesheet.csv", emit: samplesheet

    script:
    template "read_curatedMetagenomicData.py"

}