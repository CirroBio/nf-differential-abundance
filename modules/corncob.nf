process run_corncob {
    container "${params.container__corncob}"

    input:
    path "counts.csv"
    path "taxonomy.csv"
    path "metadata.csv"

    output:
    path "corncob_results.csv"

    script:
    template "run_corncob.R"
}

process split_corncob {
    container "${params.container}"
    publishDir "${params.output}/corncob/", mode: 'copy', overwrite: true

    input:
        path "corncob_results.csv"

    output:
        path "*.csv"

    script:
    template "split_corncob.py"
}

workflow corncob {
    take:
        counts
        taxonomy
        metadata

    main:
        run_corncob(
            counts,
            taxonomy,
            metadata
        )

        split_corncob(
            run_corncob.out
        )

    emit:
    split_corncob.out
}