/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    INTERSECT_BGC_ME process
    Intersects antiSMASH-lite BGC region coordinates against detected ICE/IME intervals.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process INTERSECT_BGC_ME {
    tag "${meta.id}"
    label "process_single"

    conda "${moduleDir}/../environment.yml"
    container "ghcr.io/andreassag/icescreen:v1.3.3"

    input:
    tuple val(meta), path(detected_me), path(manifest), path(antismash_json)

    output:
    tuple val(meta), path("${prefix}_bgc_me_overlaps.tsv"), emit: overlaps
    tuple val("${task.process}"), val("python"), eval("python3 --version 2>&1 | sed -n '1{s/Python //; p}'"), topic: versions, emit: versions_python


    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    intersect_bgc_me.py \\
        --sample ${meta.id} \\
        --detected-me ${detected_me} \\
        --manifest ${manifest} \\
        --antismash-json ${antismash_json} \\
        --output ${prefix}_bgc_me_overlaps.tsv
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    printf 'sample\\tme_id\\tcontig\\tme_start\\tme_end\\tbgc_region_id\\tbgc_products\\tbgc_start\\tbgc_end\\toverlap_bp\\tme_overlap_frac\\tbgc_overlap_frac\\toverlapping_locus_tags\\n' > ${prefix}_bgc_me_overlaps.tsv
    """
}
