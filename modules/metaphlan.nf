include { join_abund } from "./join_abund.nf" addParams(
    tool: "metaphlan",
    metric: params.metric
)

include { join_taxonomy } from "./join_taxonomy.nf" addParams(
    tool: "metaphlan"
)

process parse {
    container "${params.container}"
    publishDir "${params.data_output}/logs/", mode: 'copy', overwrite: true, pattern: "*.log"

    input:
    tuple val(sample), path(file)

    output:
    path "${sample}.csv", emit: abund
    path "taxonomy.csv", emit: taxonomy
    path "*.log", emit: log, optional: true

    script:
    template "read_metaphlan.py"

}

workflow read {
    take:
    input_ch

    main:

    parse(input_ch)
    join_abund(parse.out.abund.toSortedList())
    join_taxonomy(parse.out.taxonomy.toSortedList())

    emit:
    abund = join_abund.out.csv
    taxonomy = join_taxonomy.out.csv

}