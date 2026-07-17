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
    container "ghcr.io/exterex/icescreen:v1.3.3"

    input:
    tuple val(meta), path(final_sp), path(genbank)
    path db

    output:
    tuple val(meta), path("${meta.id}_detected_ME.tsv"), emit: detected_me
    tuple val(meta), path("${meta.id}_detected_SP_withMEIds.tsv"), emit: sp_with_me_ids
    tuple val(meta), path("${meta.id}_detected_ME.summary"), emit: summary
    tuple val(meta), path("${meta.id}_detected_ME.log"), emit: log
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def mode_file = "${db}/mode/${params.phylum}.yml"
    def me_conf = "${db}/icescreen.conf"
    """
    icescreen_OO \\
        -i ${final_sp} \\
        -c ${me_conf} \\
        -o ${meta.id}_detected_ME.tsv \\
        -m ${meta.id}_detected_SP_withMEIds.tsv \\
        -l ${meta.id}_detected_ME.log \\
        --gb_input ${genbank} \\
        --taxo_mode_file ${mode_file}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version 2>&1 | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}_detected_ME.tsv
    touch ${meta.id}_detected_SP_withMEIds.tsv
    touch ${meta.id}_detected_ME.summary
    touch ${meta.id}_detected_ME.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: 3.11.0
    END_VERSIONS
    """
}
