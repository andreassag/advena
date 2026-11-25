/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    DETECT_ATT_SITES process
    Detects attL/attR direct repeats and tRNA/housekeeping targets flanking integrases.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process DETECT_ATT_SITES {
    tag "${meta.id}"
    label "process_single"

    conda "${moduleDir}/../environment.yml"
    container "ghcr.io/andreassag/icescreen:v1.3.3"

    input:
    tuple val(meta), path(detected_me), path(annotated_gb), path(annotated_gff)

    output:
    tuple val(meta), path("${prefix}_att_sites.tsv"), emit: att_sites
    tuple val("${task.process}"), val("blastn"), eval("blastn -version 2>&1 | sed -n '1{s/^.*blastn: //; s/ .*//; p}' || echo 'unknown'"), topic: versions, emit: versions_blastn
    tuple val("${task.process}"), val("python"), eval("python3 --version 2>&1 | sed -n '1{s/Python //; p}'"), topic: versions, emit: versions_python


    script:
    def args = task.ext.args ?: ""
    prefix = task.ext.prefix ?: "${meta.id}"
    def gff_arg = annotated_gff ? "--gff ${annotated_gff}" : ""
    """
    detect_att_sites.py \\
        --sample ${meta.id} \\
        --detected-me ${detected_me} \\
        --genbank ${annotated_gb} \\
        ${gff_arg} \\
        ${args} \\
        --output ${prefix}_att_sites.tsv
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    printf 'sample\\tme_id\\tcontig\\tintegrase_side\\tattl_start\\tattl_end\\tattr_start\\tattr_end\\tatt_repeat_length\\tatt_identity_pct\\tatt_repeat_seq\\tatt_target_feature\\tatt_target_distance_bp\\tdetection_method\\n' > ${prefix}_att_sites.tsv
    """
}
