/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PROCESS_HMMSCAN process
    Filters hmmscan results and selects best hits.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process PROCESS_HMMSCAN {
    tag "${meta.id}"
    label "process_low"
    conda "${moduleDir}/../environment.yml"
    container "ghcr.io/exterex/icescreen:v1.3.3"

    input:
    tuple val(meta), path(hmmscan_results)
    path db

    output:
    tuple val(meta), path("${meta.id}_hmm_SP.tsv"), emit: best_hits
    tuple val(meta), path("${meta.id}_unfiltered.tsv"), emit: hmm_unfiltered
    tuple val(meta), path("${meta.id}_filtered.tsv"), emit: hmm_filtered
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def mode_file = "${db}/mode/${params.phylum}.yml"
    """
    process_hmmscan_results \\
        -i ${hmmscan_results} \\
        -c ${mode_file} \\
        --outall ${meta.id}_unfiltered.tsv \\
        --outfiltered ${meta.id}_filtered.tsv \\
        --outbest ${meta.id}_hmm_SP.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version 2>&1 | sed 's/Python //')
        pandas: \$(python3 -c "import pandas; print(pandas.__version__)")
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}_hmm_SP.tsv
    touch ${meta.id}_unfiltered.tsv
    touch ${meta.id}_filtered.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: 3.11.0
        pandas: 2.2.0
    END_VERSIONS
    """
}
