/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COHORT_SUMMARY process
    Builds cohort ME cargo summary, pan-ME eggNOG OG matrix, and MultiQC custom panels.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process COHORT_SUMMARY {
    tag "cohort_summary"
    label "process_single"

    conda "${moduleDir}/../environment.yml"
    container "ghcr.io/andreassag/icescreen:v1.3.3"

    input:
    path detected_mes, stageAs: "detected_mes/*"
    path contextualized_mes, stageAs: "context_mes/*"
    path cargo_profiles, stageAs: "cargo/*"
    path qc_summary
    path clusters_tsv

    output:
    path "cohort_me_cargo_summary.tsv", emit: cohort_summary
    path "og_presence_absence.tsv", emit: og_matrix
    path "icescreen_summary_mqc.yaml", emit: mqc_icescreen
    path "amr_card_heatmap_mqc.tsv", emit: mqc_amr
    path "cog_distribution_mqc.tsv", emit: mqc_cog
    path "mge_context_mqc.tsv", emit: mqc_mge
    tuple val("${task.process}"), val("python"), eval("python3 --version 2>&1 | sed -n '1{s/Python //; p}'"), topic: versions, emit: versions_python


    script:
    def me_arg = detected_mes ? "--detected-me ${detected_mes}" : ""
    def ctx_arg = contextualized_mes ? "--contextualized-me ${contextualized_mes}" : ""
    def cargo_arg = cargo_profiles ? "--cargo-profiles ${cargo_profiles}" : ""
    def qc_arg = qc_summary ? "--qc-summary ${qc_summary}" : ""
    def cl_arg = clusters_tsv ? "--clusters ${clusters_tsv}" : ""
    """
    build_cohort_summary.py \\
        ${me_arg} \\
        ${ctx_arg} \\
        ${cargo_arg} \\
        ${qc_arg} \\
        ${cl_arg} \\
        --outdir .
    """

    stub:
    """
    printf 'sample\\tme_id\\tcontig\\tcategory\\n' > cohort_me_cargo_summary.tsv
    printf 'og_id\\tog_level\\tcog_category\\tconsensus_description\\tn_elements\\tprevalence_pct\\tclassification\\n' > og_presence_absence.tsv
    printf "id: 'advena_icescreen_summary'\\nsection_name: 'advena: ICE/IME Summary'\\nplot_type: 'bargraph'\\ndata: {}\\n" > icescreen_summary_mqc.yaml
    printf '# id: advena_amr_card\\nSample\\tNo_AMR_Hits\\n' > amr_card_heatmap_mqc.tsv
    printf '# id: advena_cog_distribution\\nSample\\tUnassigned\\n' > cog_distribution_mqc.tsv
    printf '# id: advena_mge_context\\nSample\\tPlasmid_borne_MEs\\n' > mge_context_mqc.tsv
    """
}
