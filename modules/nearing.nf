include { filter } from "./filter.nf"

process format_nearing {
    container "${params.container}"

    input:
    tuple path("inputs/counts.csv"), path("inputs/taxonomy.csv"), path("inputs/metadata.csv")

    output:
    tuple path("counts.tsv"), path("metadata.tsv")

    script:
    template "format_nearing.py"
}

process run_nearing {
    container "${params.container__nearing}"
    publishDir "${params.data_output}/${params.method}/", mode: 'copy', overwrite: true, pattern: "*.csv"
    publishDir "${params.data_output}/logs/", mode: 'copy', overwrite: true, pattern: "*.log"

    input:
    tuple path("counts.tsv"), path("metadata.tsv")
    path script

    output:
    path "*.results.csv", emit: csv
    path "*.log", emit: log, optional: true

    script:
    template "run_nearing.sh"
}

workflow nearing {
    take:
        counts
        taxonomy
        metadata

    main:

        // Get the script that will be run
        nearing_script = file(
            "$projectDir/assets/Comparison_of_DA_microbiome_methods-1.0.2/Pipeline_scripts/Tool_scripts/Run_${params.method}.R"
        )
        if(!nearing_script.exists()){
            log.info"""
            Available methods from Nearing, et al. are:
            - Aldex2
            - ANCOM
            - Corncob
            - DESeq2
            - edgeR
            - Limma_Voom_TMM
            - Limma_Voom_TMMwsp
            - Masslin2
            - metagenomeSeq
            - t_test_rare
            - Wilcox_CLR
            - Wilcox_rare
            """
            error "Did not recognize method: ${params.method}"
        }

        filter(
            counts,
            taxonomy,
            metadata
        )

        format_nearing(filter.out.csv)

        run_nearing(format_nearing.out, nearing_script)

    emit:
        run_nearing.out.csv
}
