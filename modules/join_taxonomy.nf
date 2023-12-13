process join_taxonomy {
    container "${params.container}"
    publishDir "${params.data_output}/${params.tool}/", mode: 'copy', overwrite: true, pattern: "taxonomy.csv"
    publishDir "${params.data_output}/logs/", mode: 'copy', overwrite: true, pattern: "*.log"

    input:
    path "inputs/taxonomy.*.csv"

    output:
    path "taxonomy.csv", emit: csv
    path "*.log", emit: log

    script:
    template "join_taxonomy.py"
}