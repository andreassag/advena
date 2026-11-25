/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    DEREP_AMR_CARD process
    Dereplicates AMRFinderPlus (--plus) and CARD RGI calls by locus_tag or >=80% overlap.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process DEREP_AMR_CARD {
    tag "${meta.id}"
    label "process_single"

    conda "${moduleDir}/../environment.yml"
    container "ghcr.io/andreassag/icescreen:v1.3.3"

    input:
    tuple val(meta), path(manifest), path(me_bed), path(amrfinder_tsv), path(card_rgi_txt)

    output:
    tuple val(meta), path("${prefix}_amr_card_derep.tsv"), emit: derep_tsv
    tuple val("${task.process}"), val("python"), eval("python3 --version 2>&1 | sed -n '1{s/Python //; p}'"), topic: versions, emit: versions_python


    script:
    def args = task.ext.args ?: ""
    prefix = task.ext.prefix ?: "${meta.id}"
    def bed_arg = me_bed ? "--bed ${me_bed}" : ""
    def amr_arg = amrfinder_tsv ? "--amrfinder ${amrfinder_tsv}" : ""
    def rgi_arg = card_rgi_txt ? "--card-rgi ${card_rgi_txt}" : ""
    """
    derep_amr_card.py \\
        --sample ${meta.id} \\
        --manifest ${manifest} \\
        ${bed_arg} \\
        ${amr_arg} \\
        ${rgi_arg} \\
        ${args} \\
        --output ${prefix}_amr_card_derep.tsv
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    printf 'sample\\tme_id\\tcontig\\tlocus_tag\\tstart\\tend\\tstrand\\tcanonical_symbol\\telement_type\\telement_subtype\\tcard_aro_accession\\tcard_aro_name\\tresistance_mechanism\\tconsensus_drug_class\\tdetected_by\\trepresentative_source\\trepresentative_tier\\tamrfinder_symbol\\tamrfinder_method\\tamrfinder_identity\\tamrfinder_coverage\\tcard_cut_off\\tcard_identity\\tcard_coverage\\n' > ${prefix}_amr_card_derep.tsv
    """
}
