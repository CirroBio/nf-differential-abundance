include { join_abund as join_counts } from "./join_abund.nf" addParams(
    tool: "sourmash",
    metric: "counts"
)

include { join_abund as join_proportions } from "./join_abund.nf" addParams(
    tool: "sourmash",
    metric: "proportions"
)

include { join_taxonomy } from "./join_taxonomy.nf" addParams(
    tool: "sourmash"
)

process parse {
    container "${params.container}"
    publishDir "${params.data_output}/logs/", mode: 'copy', overwrite: true, pattern: "*.log"

    input:
    tuple val(sample), path(file)

    output:
    path "${sample}.counts.csv", emit: counts
    path "${sample}.proportions.csv", emit: proportions
    path "taxonomy.csv", emit: taxonomy
    path "*.log", emit: log, optional: true

    script:
    template "read_sourmash.py"

}

workflow read {
    take:
    input_ch

    main:

    parse(input_ch)
    join_counts(parse.out.counts.toSortedList())
    join_proportions(parse.out.proportions.toSortedList())
    join_taxonomy(parse.out.taxonomy.toSortedList())

    emit:
    counts = join_counts.out.csv
    proportions = join_proportions.out.csv
    taxonomy = join_taxonomy.out
}