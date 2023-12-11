process wilcoxon {
    container "${params.container}"
    publishDir "${params.data_output}/wilcoxon/", mode: 'copy', overwrite: true

    input:
    path "counts.csv"
    path "metadata.csv"

    output:
    path "*.csv"

    script:
    template "wilcoxon.py"
}
