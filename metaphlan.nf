#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Use the same process to parse two different metrics
// from the same set of output files
include { read as read_counts } from './modules/metaphlan.nf' addParams(metric: 'estimated_number_of_reads_from_the_clade')
include { read as read_prop } from './modules/metaphlan.nf' addParams(metric: 'relative_abundance')

// Shared workflow for running differential abundance analysis
include { differential_abundance } from './modules/differential_abundance.nf' addParams(tool: 'metaphlan')

workflow {
    if(!params.samplesheet){
        error "Must provide param: samplesheet"
    }
    if(!params.formula){
        error "Must provide param: formula"
    }

    samplesheet = file(params.samplesheet, checkIfExists: true)

    Channel
        .from(
            samplesheet
        )
        .splitCsv(
            header: true,
            sep: ","
        )
        .map {
            it -> [
                it.sample,
                file(it.file, checkIfExists: true)
            ]
        }
        .set { input_ch }

    // Parse the read counts and the proportional abundances
    read_counts(input_ch)
    read_prop(input_ch)

    // Run the differential abundance analysis
    differential_abundance(
        read_counts.out.abund,
        read_prop.out.abund,
        read_counts.out.taxonomy,
        samplesheet
    )

}