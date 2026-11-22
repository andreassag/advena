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
    container "ghcr.io/andreassag/icescreen:v1.3.3"

    input:
    tuple val(meta), path(hmmscan_results)
    path mode_file

    output:
    tuple val(meta), path("${prefix}_hmm_SP.tsv"), emit: best_hits
    tuple val(meta), path("${prefix}_unfiltered.tsv"), emit: hmm_unfiltered
    tuple val(meta), path("${prefix}_filtered.tsv"), emit: hmm_filtered
    tuple val("${task.process}"), val("icescreen"), eval("python3 -c 'import importlib.metadata; print(importlib.metadata.version(\"icescreen\"))' 2>/dev/null || echo 'unknown'"), topic: versions, emit: versions_icescreen
    tuple val("${task.process}"), val("pandas"), eval("python3 -c 'import pandas; print(pandas.__version__)' 2>/dev/null || echo 'unknown'"), topic: versions, emit: versions_pandas

    script:
    def args = task.ext.args ?: ""
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    process_hmmscan_results \\
        ${args} \\
        -i ${hmmscan_results} \\
        -c ${mode_file} \\
        --outall ${prefix}_unfiltered.tsv \\
        --outfiltered ${prefix}_filtered.tsv \\
        --outbest ${prefix}_hmm_SP.tsv
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_hmm_SP.tsv
    touch ${prefix}_unfiltered.tsv
    touch ${prefix}_filtered.tsv
    """
}
