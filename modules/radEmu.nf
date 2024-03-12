process run_radEmu {
    container "${params.container__radEmu}"
    publishDir "${params.data_output}/logs/", mode: 'copy', overwrite: true, pattern: "*.log"

    input:
    path "counts.csv"
    path "taxonomy.csv"
    path "metadata.csv"

    output:
    path "radEmu_results.csv", emit: csv
    path "*.log", emit: log, optional: true

    script:
    template "run_radEmu.R"
}

process split_radEmu {
    container "${params.container}"
    publishDir "${params.data_output}/radEmu/", mode: 'copy', overwrite: true, pattern: "*.csv"
    publishDir "${params.data_output}/logs/", mode: 'copy', overwrite: true, pattern: "*.log"

    input:
        path "radEmu_results.csv"

    output:
        path "*.csv", emit: csv
        path "*.log", emit: log, optional: true

    script:
    template "split_radEmu.py"
}

workflow radEmu {
    take:
        counts
        taxonomy
        metadata

    main:
        run_radEmu(
            counts,
            taxonomy,
            metadata
        )

        split_radEmu(
            run_radEmu.out.csv
        )

    emit:
        split_radEmu.out.csv
}
