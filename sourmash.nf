#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Read sourmash data
include { read } from './modules/sourmash.nf'

// Shared workflow for running differential abundance analysis
include { differential_abundance } from './modules/differential_abundance.nf' addParams(tool: 'sourmash')

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

    // Parse the pseudo read counts (total_weighted_hashes) and the proportional abundances
    read(input_ch)

    // Run the differential abundance analysis
    differential_abundance(
        read.out.counts,
        read.out.proportions,
        read.out.taxonomy,
        samplesheet
    )

}