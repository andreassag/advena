/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PROCESS_HMMSCAN_FP process
    Filters false positive signature proteins from merged SP hits using FP hmmscan results.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process PROCESS_HMMSCAN_FP {
    tag "${meta.id}"
    label "process_low"

    conda "${moduleDir}/../environment.yml"
    container "ghcr.io/andreassag/icescreen:v1.3.3"

    input:
    tuple val(meta), path(merged_sp), path(fp_all)
    path mode_file

    output:
    tuple val(meta), path("${prefix}_detected_SP_hmm_cleaned.tsv"), emit: sp_clean
    tuple val(meta), path("${prefix}_hmm_fp_hits.tsv"), emit: fp_hits
    tuple val("${task.process}"), val("icescreen"), eval("python3 -c 'import importlib.metadata; print(importlib.metadata.version(\"icescreen\"))' 2>/dev/null || echo 'unknown'"), topic: versions, emit: versions_icescreen
    tuple val("${task.process}"), val("pandas"), eval("python3 -c 'import pandas; print(pandas.__version__)' 2>/dev/null || echo 'unknown'"), topic: versions, emit: versions_pandas

    script:
    def args = task.ext.args ?: ""
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    process_hmmscan_fp \\
        ${args} \\
        --insp ${merged_sp} \\
        --infp ${fp_all} \\
        -c ${mode_file} \\
        --outfp ${prefix}_hmm_fp_hits.tsv \\
        --outfiltered ${prefix}_detected_SP_hmm_cleaned.tsv
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_detected_SP_hmm_cleaned.tsv
    touch ${prefix}_hmm_fp_hits.tsv
    """
}
