#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Use the same process to parse two different metrics
// from the same set of output files
include { read as read_counts } from './modules/metaphlan.nf' addParams(metric: 'estimated_number_of_reads_from_the_clade')
include { read as read_prop } from './modules/metaphlan.nf' addParams(metric: 'relative_abundance')
include { ad_metaphlan } from './modules/ad_metaphlan.nf'

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

    // Parse the read counts and the proportional abundances
    read_counts(input_ch)
    read_prop(input_ch)

    if(params.method == "corncob"){

        // Run corncob for stats
        corncob(
            read_counts.out.abund,
            read_counts.out.taxonomy,
            samplesheet
        )
        stats_output = corncob.out

    } else {
        if(params.method != "wilcoxon"){
            error "Parameter 'method' not recognized: ${params.method}"
        }
        // Run wilcoxon for stats
        wilcoxon(
            read_counts.out.abund,
            samplesheet
        )
        stats_output = wilcoxon.out

    }

    // Make an AnnData object with both proportions and counts
    // as well as the stats results
    ad_metaphlan(
        read_counts.out.abund,
        read_prop.out.abund,
        read_prop.out.taxonomy,
        samplesheet,
        stats_output
    )

    // Make the visualization elements
    viz(ad_metaphlan.out)
}