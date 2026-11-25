/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    EGGNOG_DOWNLOAD process
    Downloads and caches the eggNOG-mapper v5 reference database (~50 GB) when
    --auto_download_large_dbs true is explicitly enabled.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process EGGNOG_DOWNLOAD {
    tag "eggnog_db"
    label "process_low"

    conda "bioconda::eggnog-mapper=2.1.12"
    container "quay.io/biocontainers/eggnog-mapper:2.1.12--pyhdfd78af_0"

    output:
    path "eggnog_db", emit: db
    tuple val("${task.process}"), val("eggnogmapper"), eval("emapper.py --version 2>&1 | sed -n '1{s/emapper-//;p}' || echo 'unknown'"), topic: versions, emit: versions_eggnogmapper


    script:
    """
    mkdir -p eggnog_db
    download_eggnog_data.py -y --data_dir eggnog_db
    """

    stub:
    """
    mkdir -p eggnog_db
    touch eggnog_db/eggnog.db eggnog_db/eggnog.taxa.db eggnog_db/eggnog_proteins.dmnd
    """
}
