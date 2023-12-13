#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Read a relative abundance table
include { read } from './modules/table.nf'

// Shared workflow for running differential abundance analysis
include { differential_abundance } from './modules/differential_abundance.nf' addParams(tool: 'table')

workflow {
    if(!params.input){
        error "Must provide param: input"
    }
    if(!params.metadata){
        error "Must provide param: metadata"
    }
    if(!params.formula){
        error "Must provide param: formula"
    }

    input = file(params.input, checkIfExists: true)
    metadata = file(params.metadata, checkIfExists: true)

    // Parse the abundances and the sample metadata
    read(
        input,
        metadata
    )

    // Run the differential abundance analysis
    differential_abundance(
        read.out.counts,
        read.out.proportions,
        read.out.taxonomy,
        read.out.samplesheet
    )

}