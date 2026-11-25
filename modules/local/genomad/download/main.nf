/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    GENOMAD_DOWNLOAD process
    Downloads and caches the geNomad plasmid and virus/prophage reference database.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process GENOMAD_DOWNLOAD {
    tag "genomad_db"
    label "process_low"

    conda "bioconda::genomad=1.8.1"
    container "quay.io/biocontainers/genomad:1.8.1--pyhdfd78af_0"

    output:
    path "genomad_db", emit: db
    tuple val("${task.process}"), val("genomad"), eval("genomad --version 2>&1 | sed -n '1{s/geNomad, version //; p}' || echo 'unknown'"), topic: versions, emit: versions_genomad


    script:
    """
    genomad download-database .
    """

    stub:
    """
    mkdir -p genomad_db
    touch genomad_db/version.txt
    """
}
