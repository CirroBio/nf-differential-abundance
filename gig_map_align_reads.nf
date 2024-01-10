#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Read output from gig-map/align_reads.nf
include { read } from './modules/gig_map_align_reads.nf'

// Shared workflow for running differential abundance analysis
include { differential_abundance } from './modules/differential_abundance.nf' addParams(tool: 'gig_map_align_reads')

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

    read_alignments = file(params.input, checkIfExists: true)
    metadata = file(params.metadata, checkIfExists: true)

    // Parse the abundances and the sample metadata
    read(
        read_alignments,
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