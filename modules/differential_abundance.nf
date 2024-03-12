// Functions shared across input types
include { make_anndata } from './make_anndata.nf'
include { corncob } from './corncob.nf'
include { radEmu } from './radEmu.nf'
include { mannwhitneyu } from './mannwhitneyu.nf'
include { viz } from './viz.nf'

workflow differential_abundance {

    take:
        counts
        proportions
        taxonomy
        samplesheet

    main:

        if (params.method == "corncob") {

            // Run corncob for stats
            corncob(
                counts,
                taxonomy,
                samplesheet
            )
            stats_output = corncob.out

        } else if (params.method == "radEmu") {
            // Run radEmu for stats
            radEmu(
                counts,
                taxonomy,
                samplesheet
            )
            stats_output = radEmu.out

        } else if (params.method == "mannwhitneyu") {
            // Run mannwhitneyu for stats
            mannwhitneyu(
                proportions,
                samplesheet
            )
            stats_output = mannwhitneyu.out.csv

        } else {
            error "Parameter 'method' not recognized: ${params.method}"
        }

        // Make an AnnData object with both proportions and counts
        // as well as the stats results
        make_anndata(
            counts,
            proportions,
            taxonomy,
            samplesheet,
            stats_output
        )

        // Make the visualization elements
        viz(make_anndata.out[0])

}