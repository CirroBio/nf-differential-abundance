process vitessce {
    container "${params.container}"
    publishDir "${params.web_output}", mode: 'copy', overwrite: true

    input:
        tuple path("adata.h5ad"), path("input_config.json")

    output:
        path "*.zarr", hidden: true
        path "*.json"

    script:
    template "ad_vitessce.py"

}

process normalize {
    container "${params.container}"
    publishDir "${params.data_output}/logs/", mode: 'copy', overwrite: true, pattern: "*.log"

    input:
        tuple path("adata.h5ad"), path("input_config.json")

    output:
        tuple path("normalized.h5ad"), path("output_config.json")
        path "*.log", optional: true

    script:
    template "ad_normalize.py"
}

process sort {
    container "${params.container}"
    publishDir "${params.data_output}/logs/", mode: 'copy', overwrite: true, pattern: "*.log"

    input:
        tuple path("adata.h5ad"), path("input_config.json")

    output:
        tuple path("sorted.h5ad"), path("output_config.json")
        path "*.log", optional: true

    script:
    template "ad_sort.py"
}

process annotate {
    container "${params.container}"
    publishDir "${params.data_output}/logs/", mode: 'copy', overwrite: true, pattern: "*.log"
    publishDir "${params.data_output}/anndata/", mode: 'copy', overwrite: true, pattern: "*.h5ad"
    publishDir "${params.data_output}/csv/", mode: 'copy', overwrite: true, pattern: "*.csv"

    input:
        tuple path("adata.h5ad"), path("input_config.json")

    output:
        tuple path("annotated.h5ad"), path("output_config.json")
        path "*.log", optional: true
        path "*.csv"

    script:
    template "ad_annotate.py"
}

workflow viz {
    take:
        anndata

    main:
        normalize(anndata)
        sort(normalize.out[0])
        annotate(sort.out[0])
        vitessce(annotate.out[0])
}