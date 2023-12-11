#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Read output from curatedMetagenomicDataTerminal
include { read } from './modules/curatedMetagenomicData.nf'

// Shared workflow for running differential abundance analysis
include { differential_abundance } from './modules/differential_abundance.nf' addParams(tool: 'curatedMetagenomicData')

workflow {
    if(!params.input){
        error "Must provide param: input"
    }
    if(!params.formula){
        error "Must provide param: formula"
    }

    input = file(params.input, checkIfExists: true)

    // Parse the proportional abundances and the sample metadata
    read(input)

    // Run the differential abundance analysis
    differential_abundance(
        read.out.counts,
        read.out.proportions,
        read.out.taxonomy,
        read.out.samplesheet
    )

}