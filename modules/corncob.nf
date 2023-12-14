process run_corncob {
    container "${params.container__corncob}"
    publishDir "${params.data_output}/logs/", mode: 'copy', overwrite: true, pattern: "*.log"

    input:
    path "counts.csv"
    path "taxonomy.csv"
    path "metadata.csv"

    output:
    path "corncob_results.csv", emit: csv
    path "*.log", emit: log, optional: true

    script:
    template "run_corncob.R"
}

process split_corncob {
    container "${params.container}"
    publishDir "${params.data_output}/corncob/", mode: 'copy', overwrite: true, pattern: "*.csv"
    publishDir "${params.data_output}/logs/", mode: 'copy', overwrite: true, pattern: "*.log"

    input:
        path "corncob_results.csv"

    output:
        path "*.csv", emit: csv
        path "*.log", emit: log, optional: true

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
            run_corncob.out.csv
        )

    emit:
        split_corncob.out.csv
}
