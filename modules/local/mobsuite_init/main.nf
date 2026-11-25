/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    MOBSUITE_INIT process
    Downloads and caches the MOB-suite plasmid reference database.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process MOBSUITE_INIT {
    tag "mobsuite_db"
    label "process_single"

    conda "bioconda::mob_suite=3.1.9"
    container "quay.io/biocontainers/mob_suite:3.1.9--pyhdfd78af_1"

    output:
    path "mobsuite_db", emit: db
    tuple val("${task.process}"), val("mobsuite"), eval("mob_init --version 2>&1 | sed -n '1{s/mob_init //; p}' || echo 'unknown'"), topic: versions, emit: versions_mobsuite


    script:
    """
    mkdir -p mobsuite_db
    mob_init --database_directory mobsuite_db
    """

    stub:
    """
    mkdir -p mobsuite_db
    touch mobsuite_db/status.txt
    """
}
