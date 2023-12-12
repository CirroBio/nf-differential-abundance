process mannwhitneyu {
    container "${params.container}"
    publishDir "${params.data_output}/mannwhitneyu/", mode: 'copy', overwrite: true

    input:
    path "counts.csv"
    path "metadata.csv"

    output:
    path "*.csv"

    script:
    template "mannwhitneyu.py"
}
