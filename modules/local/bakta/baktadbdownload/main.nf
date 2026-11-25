/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    BAKTA_BAKTADBDOWNLOAD process
    Downloads and caches the Bakta reference database (light or full).
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process BAKTA_BAKTADBDOWNLOAD {
    tag "${db_type}"
    label "process_low"

    conda "bioconda::bakta=1.10.4"
    container "quay.io/biocontainers/bakta:1.10.4--pyhdfd78af_0"

    input:
    val db_type

    output:
    path "bakta_db", emit: db
    tuple val("${task.process}"), val("bakta"), eval("bakta --version 2>&1 | sed -n '1{s/bakta //; p}' || echo 'unknown'"), topic: versions, emit: versions_bakta


    script:
    def args = task.ext.args ?: ""
    """
    bakta_db download \\
        --output bakta_db \\
        --type ${db_type} \\
        ${args}
    """

    stub:
    """
    mkdir -p bakta_db/db-${db_type}
    touch bakta_db/db-${db_type}/bakta.db
    """
}
