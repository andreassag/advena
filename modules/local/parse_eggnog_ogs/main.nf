/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PARSE_EGGNOG_OGS process
    Extracts Bacteria (@2) and Bacillota (@1239) OGs, COG categories, and cargo roles.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process PARSE_EGGNOG_OGS {
    tag "${meta.id}"
    label "process_single"

    conda "${moduleDir}/../environment.yml"
    container "ghcr.io/andreassag/icescreen:v1.3.3"

    input:
    tuple val(meta), path(emapper_annotations), path(manifest)

    output:
    tuple val(meta), path("${prefix}_cargo_ogs.tsv"), emit: cargo_ogs
    tuple val("${task.process}"), val("python"), eval("python3 --version 2>&1 | sed -n '1{s/Python //; p}'"), topic: versions, emit: versions_python


    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    parse_eggnog_ogs.py \\
        --sample ${meta.id} \\
        --emapper ${emapper_annotations} \\
        --manifest ${manifest} \\
        --output ${prefix}_cargo_ogs.tsv
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    printf 'sample\\tme_id\\tcontig\\tlocus_tag\\tis_signature_protein\\tcargo_role\\troot_cog\\tbacteria_og\\tbacillota_og\\tnarrowest_og\\tcog_category\\tpreferred_name\\tdescription\\tkegg_ko\\tkegg_pathway\\tpfams\\tgos\\tseed_ortholog\\tevalue\\tscore\\n' > ${prefix}_cargo_ogs.tsv
    """
}
