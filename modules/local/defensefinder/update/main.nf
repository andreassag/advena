/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    DEFENSEFINDER_UPDATE process
    Downloads and caches DefenseFinder (MacSyFinder + CasFinder) models.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process DEFENSEFINDER_UPDATE {
    tag "defensefinder_db"
    label "process_single"

    conda "bioconda::defense-finder=1.3.0"
    container "quay.io/biocontainers/defense-finder:1.3.0--pyhdfd78af_0"

    output:
    path "defensefinder_models", emit: db
    tuple val("${task.process}"), val("defensefinder"), eval("defense-finder version 2>&1 | sed -n '1p' || echo 'unknown'"), topic: versions, emit: versions_defensefinder


    script:
    """
    mkdir -p defensefinder_models
    defense-finder update --models-dir defensefinder_models
    """

    stub:
    """
    mkdir -p defensefinder_models/defense-finder-models
    touch defensefinder_models/defense-finder-models/metadata.yml
    """
}
