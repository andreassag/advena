/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONTEXTUALIZE_ME process
    Overlaps MOB-suite, geNomad, ISEScan, IntegronFinder, and attL/attR calls with ICE/IMEs.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process CONTEXTUALIZE_ME {
    tag "${meta.id}"
    label "process_single"

    conda "${moduleDir}/../environment.yml"
    container "ghcr.io/andreassag/icescreen:v1.3.3"

    input:
    tuple val(meta), path(detected_me), path(mobsuite_report), path(genomad_plasmid), path(genomad_virus), path(isescan_tsv), path(integron_file), path(att_sites_tsv)

    output:
    tuple val(meta), path("${prefix}_detected_ME_contextualized.tsv"), emit: contextualized_me
    tuple val("${task.process}"), val("python"), eval("python3 --version 2>&1 | sed -n '1{s/Python //; p}'"), topic: versions, emit: versions_python


    script:
    def args = task.ext.args ?: ""
    prefix = task.ext.prefix ?: "${meta.id}"
    def mob_arg = mobsuite_report ? "--mobsuite ${mobsuite_report}" : ""
    def gp_arg = genomad_plasmid ? "--genomad-plasmid ${genomad_plasmid}" : ""
    def gv_arg = genomad_virus ? "--genomad-virus ${genomad_virus}" : ""
    def is_arg = isescan_tsv ? "--isescan ${isescan_tsv}" : ""
    def int_arg = integron_file ? "--integronfinder ${integron_file}" : ""
    def att_arg = att_sites_tsv ? "--att-sites ${att_sites_tsv}" : ""
    """
    contextualize_me.py \\
        --detected-me ${detected_me} \\
        ${mob_arg} \\
        ${gp_arg} \\
        ${gv_arg} \\
        ${is_arg} \\
        ${int_arg} \\
        ${att_arg} \\
        ${args} \\
        --output ${prefix}_detected_ME_contextualized.tsv
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    printf 'ICE_IME_id\\tGenome_accession\\tCategory_of_element\\tStart_of_most_upstream_SP\\tStop_of_most_downstream_SP\\tplasmid_borne_warning\\tmobsuite_molecule_type\\tmobsuite_cluster_id\\tgenomad_plasmid_score\\tcomposite_prophage_ids\\tflanking_is_families\\tinternal_is_count\\tintegron_ids\\tatt_target_feature\\tatt_repeat_length\\tatt_repeat_seq\\n' > ${prefix}_detected_ME_contextualized.tsv
    """
}
