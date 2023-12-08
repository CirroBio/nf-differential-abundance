process parse {
    container "${params.container}"

    input:
    tuple val(sample), path(file)

    output:
    path "${sample}.csv", emit: abund
    path "taxonomy.csv", emit: taxonomy

    script:
    template "read_metaphlan.py"

}

process join_abund {
    container "${params.container}"
    publishDir "${params.output}/metaphlan/",
        mode: 'copy',
        overwrite: true,
        saveAs: { filename -> "${params.tax_level}.${params.metric}.csv" }

    input:
    path "inputs/"

    output:
    path "abund.csv"

    script:
    template "join_abund.py"

}

process join_taxonomy {
    container "${params.container}"
    publishDir "${params.output}/metaphlan/", mode: 'copy', overwrite: true

    input:
    path "inputs/taxonomy.*.csv"

    output:
    path "taxonomy.csv"

    script:
    template "join_taxonomy.py"

}

workflow read {
    take:
    input_ch

    main:

    parse(input_ch)
    join_abund(parse.out.abund.toSortedList())
    join_taxonomy(parse.out.taxonomy.toSortedList())

    emit:
    abund = join_abund.out
    taxonomy = join_taxonomy.out

}