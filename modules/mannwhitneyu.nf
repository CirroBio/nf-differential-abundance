process mannwhitneyu {
    container "${params.container}"
    publishDir "${params.data_output}/mannwhitneyu/", mode: 'copy', overwrite: true, pattern: "*.csv"
    publishDir "${params.data_output}/logs/", mode: 'copy', overwrite: true, pattern: "*.log"

    input:
    path "counts.csv"
    path "metadata.csv"

    output:
    path "*.csv", emit: csv
    path "*.log", emit: log, optional: true

    script:
    template "mannwhitneyu.py"
}
