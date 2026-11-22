/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    REANNOT_SP process
    Re-annotates signature proteins (XerS classification).
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process REANNOT_SP {
    tag "${meta.id}"
    label "process_low"

    conda "${moduleDir}/../environment.yml"
    container "ghcr.io/andreassag/icescreen:v1.3.3"

    input:
    tuple val(meta), path(sp_clean), path(intyr_hits)
    path mode_file

    output:
    tuple val(meta), path("${prefix}_detected_SP.tsv"), emit: final_sp
    tuple val(meta), path("${prefix}_detected_SP_hmm_cleaned_reannotated.tsv"), emit: reannotated
    tuple val("${task.process}"), val("icescreen"), eval("python3 -c 'import importlib.metadata; print(importlib.metadata.version(\"icescreen\"))' 2>/dev/null || echo 'unknown'"), topic: versions, emit: versions_icescreen
    tuple val("${task.process}"), val("pandas"), eval("python3 -c 'import pandas; print(pandas.__version__)' 2>/dev/null || echo 'unknown'"), topic: versions, emit: versions_pandas

    script:
    def args = task.ext.args ?: ""
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    reannot_SP \\
        ${args} \\
        -a ${sp_clean} \\
        -b ${intyr_hits} \\
        -c ${mode_file} \\
        --outall ${prefix}_detected_SP_hmm_cleaned_reannotated.tsv \\
        --outfiltered ${prefix}_detected_SP.tsv
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_detected_SP.tsv
    touch ${prefix}_detected_SP_hmm_cleaned_reannotated.tsv
    """
}
