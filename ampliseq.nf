#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Read output from nf-core/ampliseq
include { read } from './modules/ampliseq.nf'

// Shared workflow for running differential abundance analysis
include { differential_abundance } from './modules/differential_abundance.nf' addParams(tool: 'ampliseq')

workflow {
    if(!params.counts){
        error "Must provide param: counts"
    }
    if(!params.taxonomy){
        error "Must provide param: taxonomy"
    }
    if(!params.metadata){
        error "Must provide param: metadata"
    }
    if(!params.formula){
        error "Must provide param: formula"
    }

    counts = file(params.counts, checkIfExists: true)
    taxonomy = file(params.taxonomy, checkIfExists: true)
    metadata = file(params.metadata, checkIfExists: true)

    // Parse the abundances and the sample metadata
    read(
        counts,
        taxonomy,
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