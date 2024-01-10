process read {
    container "${params.container}"
    publishDir "${params.data_output}/logs/", mode: 'copy', overwrite: true, pattern: "*.log"

    input:
    file "read_alignments.csv.gz"
    file "metadata.csv"

    output:
    path "counts.csv", emit: counts
    path "proportions.csv", emit: proportions
    path "taxonomy.csv", emit: taxonomy
    path "samplesheet.csv", emit: samplesheet
    path "*log", emit: log

    script:
    template "read_gig_map_align_reads.py"

}