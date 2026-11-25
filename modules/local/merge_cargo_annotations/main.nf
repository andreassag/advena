/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    MERGE_CARGO_ANNOTATIONS process
    Joins AMR/CARD, eggNOG OGs, DefenseFinder, and BGC annotations into master TSV & GFF3.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process MERGE_CARGO_ANNOTATIONS {
    tag "${meta.id}"
    label "process_single"

    conda "${moduleDir}/../environment.yml"
    container "ghcr.io/andreassag/icescreen:v1.3.3"

    input:
    tuple val(meta), path(manifest), path(icescreen_gff), path(amr_tsv), path(eggnog_tsv), path(defense_genes), path(defense_systems), path(bgc_tsv)

    output:
    tuple val(meta), path("${prefix}_me_cargo_profile.tsv"), emit: cargo_profile
    tuple val(meta), path("${prefix}_advena.gff.gz"), emit: enriched_gff
    tuple val("${task.process}"), val("python"), eval("python3 --version 2>&1 | sed -n '1{s/Python //; p}'"), topic: versions, emit: versions_python


    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    def gff_arg = icescreen_gff ? "--icescreen-gff ${icescreen_gff}" : ""
    def amr_arg = amr_tsv ? "--amr ${amr_tsv}" : ""
    def egg_arg = eggnog_tsv ? "--eggnog ${eggnog_tsv}" : ""
    def def_g_arg = defense_genes ? "--defense-genes ${defense_genes}" : ""
    def def_s_arg = defense_systems ? "--defense-systems ${defense_systems}" : ""
    def bgc_arg = bgc_tsv ? "--bgc ${bgc_tsv}" : ""
    """
    merge_cargo_annotations.py \\
        --manifest ${manifest} \\
        ${gff_arg} \\
        ${amr_arg} \\
        ${egg_arg} \\
        ${def_g_arg} \\
        ${def_s_arg} \\
        ${bgc_arg} \\
        --out-profile ${prefix}_me_cargo_profile.tsv \\
        --out-gff ${prefix}_advena.gff.gz
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    printf 'sample\\tme_id\\tme_category\\tcontig\\tlocus_tag\\tprotein_id\\torig_start\\torig_end\\tstrand\\tis_signature_protein\\tsp_type\\tproduct\\tamr_gene\\tamr_element_type\\tamr_drug_class\\tcard_aro\\tamr_detected_by\\teggnog_bacteria_og\\teggnog_bacillota_og\\teggnog_narrowest_og\\tcog_category\\teggnog_description\\tcargo_role\\tdefense_system\\tdefense_gene\\tbgc_type\\n' > ${prefix}_me_cargo_profile.tsv
    printf '##gff-version 3\\n' | gzip > ${prefix}_advena.gff.gz
    """
}
