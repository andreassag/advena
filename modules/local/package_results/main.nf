/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PACKAGE_RESULTS process
    Archives intermediate files and generates param.conf.gz run summary.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process PACKAGE_RESULTS {
    tag "${meta.id}"
    label "process_low"
    conda "${moduleDir}/../environment.yml" // TODO: Might not need the conda and container environment anymore for this process
    container "ghcr.io/andreassag/icescreen:v1.3.3"

    input:
    tuple val(meta),
        path(faa),
        path(blast_raw),
        path(blast_best),
        path(blast_unfiltered),
        path(blast_filtered),
        path(blast_filtered_best),
        path(hmm_compiled),
        path(hmm_per_profile_tsv),
        path(hmm_per_profile_txt),
        path(hmm_best),
        path(hmm_unfiltered),
        path(hmm_filtered),
        path(merged_sp),
        path(sp_no_fp),
        path(sp_faa),
        path(fp_all),
        path(fp_hits),
        path(fp_per_profile_tsv),
        path(fp_per_profile_txt),
        path(sp_reannotated),
        path(final_sp),
        path(me_log)
    path db

    output:
    tuple val(meta), path("tmp_intermediate_files.tar.gz"), emit: archive
    tuple val(meta), path("param.conf.gz"), emit: param_conf
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def mode_file = "${db}/mode/${params.phylum}.yml"
    def S = "${meta.id}"
    def R = "results/${meta.id}"
    """
    # Build the intermediate file directory tree mirroring the original ICEscreen structure
    mkdir -p ${R}/icescreen_detection_SP/Blast_mode/blastp_output
    mkdir -p ${R}/icescreen_detection_SP/Blast_mode/filtered_results
    mkdir -p ${R}/icescreen_detection_SP/Blast_mode/unfiltered_results
    mkdir -p ${R}/icescreen_detection_SP/HMM_mode/hmmscan_output
    mkdir -p ${R}/icescreen_detection_SP/HMM_mode/filtered_results
    mkdir -p ${R}/icescreen_detection_SP/HMM_mode/unfiltered_results
    mkdir -p ${R}/icescreen_detection_SP/hits_cleaning/proteins_to_remove/hmmscan_compiled_output
    mkdir -p ${R}/icescreen_detection_SP/hits_cleaning/proteins_to_remove/hmmscan_output

    # Top-level: FAA file
    cp ${faa} ${R}/

    # icescreen_detection_SP/: final detected SP file
    cp ${final_sp} ${R}/icescreen_detection_SP/

    # Blast_mode/
    cp ${blast_best} ${R}/icescreen_detection_SP/Blast_mode/
    cp ${blast_raw} ${R}/icescreen_detection_SP/Blast_mode/blastp_output/
    cp ${blast_filtered} ${R}/icescreen_detection_SP/Blast_mode/filtered_results/
    cp ${blast_filtered_best} ${R}/icescreen_detection_SP/Blast_mode/filtered_results/
    cp ${blast_unfiltered} ${R}/icescreen_detection_SP/Blast_mode/unfiltered_results/

    # HMM_mode/
    cp ${hmm_best} ${R}/icescreen_detection_SP/HMM_mode/
    cp ${hmm_compiled} ${R}/icescreen_detection_SP/HMM_mode/hmmscan_output/${S}.tsv
    cp ${hmm_per_profile_tsv} ${R}/icescreen_detection_SP/HMM_mode/hmmscan_output/
    cp ${hmm_per_profile_txt} ${R}/icescreen_detection_SP/HMM_mode/hmmscan_output/
    cp ${hmm_unfiltered} ${R}/icescreen_detection_SP/HMM_mode/unfiltered_results/
    cp ${hmm_filtered} ${R}/icescreen_detection_SP/HMM_mode/filtered_results/

    # hits_cleaning/
    cp ${merged_sp} ${R}/icescreen_detection_SP/hits_cleaning/
    cp ${sp_faa} ${R}/icescreen_detection_SP/hits_cleaning/
    cp ${sp_no_fp} ${R}/icescreen_detection_SP/hits_cleaning/
    cp ${sp_reannotated} ${R}/icescreen_detection_SP/hits_cleaning/
    cp ${fp_hits} ${R}/icescreen_detection_SP/hits_cleaning/proteins_to_remove/
    cp ${fp_all} ${R}/icescreen_detection_SP/hits_cleaning/proteins_to_remove/hmmscan_compiled_output/${S}.tsv
    cp ${fp_per_profile_tsv} ${R}/icescreen_detection_SP/hits_cleaning/proteins_to_remove/hmmscan_output/
    cp ${fp_per_profile_txt} ${R}/icescreen_detection_SP/hits_cleaning/proteins_to_remove/hmmscan_output/

    tar -czf tmp_intermediate_files.tar.gz results/

    # Generate param.conf: ME detection log params + mode YAML content
    tail -n +8 ${me_log} > param.conf
    tail -n +2 ${mode_file} >> param.conf
    gzip param.conf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tar: \$(tar --version | head -1)
        gzip: \$(gzip --version | head -1)
    END_VERSIONS
    """

    stub:
    """
    touch tmp_intermediate_files.tar.gz
    echo "" | gzip > param.conf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tar: tar (GNU tar) 1.34
        gzip: gzip 1.12
    END_VERSIONS
    """
}
