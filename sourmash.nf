#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Use the same process to parse two different metrics
// from the same set of output files
include { read } from './modules/sourmash.nf'
include { make_anndata } from './modules/make_anndata.nf' addParams(tool: 'sourmash')

// Functions shared across input types
include { corncob } from './modules/corncob.nf'
include { wilcoxon } from './modules/wilcoxon.nf'
include { viz } from './modules/viz.nf'

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

    if(params.method == "corncob"){

        // Run corncob for stats
        corncob(
            read.out.counts,
            read.out.taxonomy,
            samplesheet
        )
        stats_output = corncob.out

    } else {
        if(params.method != "wilcoxon"){
            error "Parameter 'method' not recognized: ${params.method}"
        }
        // Run wilcoxon for stats
        wilcoxon(
            read.out.proportions,
            samplesheet
        )
        stats_output = wilcoxon.out

    }

    // Make an AnnData object with both proportions and counts
    // as well as the stats results
    make_anndata(
        read.out.counts,
        read.out.proportions,
        read.out.taxonomy,
        samplesheet,
        stats_output
    )

    // Make the visualization elements
    viz(make_anndata.out)
}