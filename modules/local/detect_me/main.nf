/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    DETECT_ME process
    Detects ICE and IME mobile element structures from signature proteins.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process DETECT_ME {
    tag "${meta.id}"
    label "process_low"

    conda "${moduleDir}/../environment.yml"
    container "ghcr.io/andreassag/icescreen:v1.3.3"

    input:
    tuple val(meta), path(final_sp), path(genbank)
    path mode_file
    path me_conf

    output:
    tuple val(meta), path("${prefix}_detected_ME.tsv"), emit: detected_me
    tuple val(meta), path("${prefix}_detected_SP_withMEIds.tsv"), emit: sp_with_me_ids
    tuple val(meta), path("${prefix}_detected_ME.summary"), emit: summary
    tuple val(meta), path("${prefix}_detected_ME.log"), emit: log
    tuple val("${task.process}"), val("icescreen"), eval("python3 -c 'import importlib.metadata; print(importlib.metadata.version(\"icescreen\"))' 2>/dev/null || echo 'unknown'"), topic: versions, emit: versions_icescreen
    tuple val("${task.process}"), val("python"), eval("python3 --version 2>&1 | sed 's/Python //'"), topic: versions, emit: versions_python

    script:
    def args = task.ext.args ?: ""
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    icescreen_OO \\
        ${args} \\
        -i ${final_sp} \\
        -c ${me_conf} \\
        -o ${prefix}_detected_ME.tsv \\
        -m ${prefix}_detected_SP_withMEIds.tsv \\
        -l ${prefix}_detected_ME.log \\
        --gb_input ${genbank} \\
        --taxo_mode_file ${mode_file}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_detected_ME.tsv
    touch ${prefix}_detected_SP_withMEIds.tsv
    touch ${prefix}_detected_ME.summary
    touch ${prefix}_detected_ME.log
    """
}
