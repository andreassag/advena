/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    AMRFINDERPLUS_UPDATE process
    Downloads and caches the NCBI AMRFinderPlus reference database.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process AMRFINDERPLUS_UPDATE {
    tag "amrfinderplus_db"
    label "process_single"

    conda "bioconda::ncbi-amrfinderplus=4.0.19"
    container "quay.io/biocontainers/ncbi-amrfinderplus:4.0.19--h6e70893_0"

    output:
    path "amrfinderplus_db", emit: db
    tuple val("${task.process}"), val("amrfinderplus"), eval("amrfinder --version 2>&1 | sed -n '1p' || echo 'unknown'"), topic: versions, emit: versions_amrfinderplus


    script:
    """
    mkdir -p amrfinderplus_db
    amrfinder_update --force_update --database amrfinderplus_db
    """

    stub:
    """
    mkdir -p amrfinderplus_db
    touch amrfinderplus_db/version.txt
    """
}
