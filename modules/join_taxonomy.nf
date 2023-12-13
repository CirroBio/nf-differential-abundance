process join_taxonomy {
    container "${params.container}"
    publishDir "${params.data_output}/${params.tool}/", mode: 'copy', overwrite: true

    input:
    path "inputs/taxonomy.*.csv"

    output:
    path "taxonomy.csv"

    script:
    template "join_taxonomy.py"
}