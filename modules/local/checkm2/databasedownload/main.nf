/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CHECKM2_DATABASEDOWNLOAD process
    Downloads and caches the CheckM2 DIAMOND reference database.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process CHECKM2_DATABASEDOWNLOAD {
    tag "checkm2_db"
    label "process_low"

    conda "bioconda::checkm2=1.0.2"
    container "quay.io/biocontainers/checkm2:1.0.2--pyh7cba7a3_0"

    output:
    path "checkm2_db", emit: db
    tuple val("${task.process}"), val("checkm2"), eval("checkm2 --version 2>&1 | sed -n '1p' || echo 'unknown'"), topic: versions, emit: versions_checkm2


    script:
    def args = task.ext.args ?: ""
    """
    mkdir -p checkm2_db
    checkm2 database --download --path checkm2_db ${args}
    """

    stub:
    """
    mkdir -p checkm2_db
    touch checkm2_db/uniref100.KO.1.dmnd
    """
}
