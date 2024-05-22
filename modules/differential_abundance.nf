// Functions shared across input types
include { make_anndata } from './make_anndata.nf'
include { corncob } from './corncob.nf'
include { radEmu } from './radEmu.nf'
include { mannwhitneyu } from './mannwhitneyu.nf'
include { nearing } from './nearing.nf'
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
            // For any other method, try to run it using the code published
            // by Nearing, et al.
            nearing(
                counts,
                taxonomy,
                samplesheet
            )
            stats_output = nearing.out
        }

        // Skip the downstream for several methods, which do not provide
        // a p-value or estimated coefficient of association
        if(!["ANCOM", "Corncob", "t_test_rare", "Wilcox_CLR", "Wilcox_rare"].contains(params.method)){

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
            if(params.run_viz){
                viz(make_anndata.out[0])
            }

        }

}