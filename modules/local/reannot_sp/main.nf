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
    container "ghcr.io/exterex/icescreen:v1.3.3"

    input:
    tuple val(meta), path(sp_clean), path(intyr_hits)
    path db

    output:
    tuple val(meta), path("${meta.id}_detected_SP.tsv"), emit: final_sp
    tuple val(meta), path("${meta.id}_detected_SP_hmm_cleaned_reannotated.tsv"), emit: reannotated
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def mode_file = "${db}/mode/${params.phylum}.yml"
    """
    reannot_SP \\
        -a ${sp_clean} \\
        -b ${intyr_hits} \\
        -c ${mode_file} \\
        --outall ${meta.id}_detected_SP_hmm_cleaned_reannotated.tsv \\
        --outfiltered ${meta.id}_detected_SP.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version 2>&1 | sed 's/Python //')
        pandas: \$(python3 -c "import pandas; print(pandas.__version__)")
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}_detected_SP.tsv
    touch ${meta.id}_detected_SP_hmm_cleaned_reannotated.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: 3.11.0
        pandas: 2.2.0
    END_VERSIONS
    """
}
