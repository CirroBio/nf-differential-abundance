process vitessce {
    container "${params.container}"
    publishDir "${params.web_output}", mode: 'copy', overwrite: true

    input:
        tuple path("adata.h5ad"), path("config.json")

    output:
        path "*.zarr", hidden: true

    script:
    template "ad_vitessce.py"

}

process normalize {
    container "${params.container}"

    input:
        tuple path("adata.h5ad"), path("config.json")

    output:
        tuple path("normalized.h5ad"), path("config.json")

    script:
    template "ad_normalize.py"
}

process sort {
    container "${params.container}"

    input:
        tuple path("adata.h5ad"), path("config.json")

    output:
        tuple path("sorted.h5ad"), path("config.json")

    script:
    template "ad_sort.py"
}

process annotate {
    container "${params.container}"

    input:
        tuple path("adata.h5ad"), path("config.json")

    output:
        tuple path("annotated.h5ad"), path("config.json")

    script:
    template "ad_annotate.py"
}

workflow viz {
    take:
        anndata

    main:
        anndata \
        | normalize \
        | sort \
        | annotate \
        | vitessce
}