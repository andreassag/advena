/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PACKAGE_RESULTS process
    Archives intermediate files and generates param.conf.gz run summary.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process PACKAGE_RESULTS {
    tag "${meta.id}"
    label "process_single"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/52/52ccce28d2ab928ab862e25aae26314d69c8e38bd41ca9431c67ef05221348aa/data'
        : 'community.wave.seqera.io/library/coreutils_grep_gzip_lbzip2_pruned:838ba80435a629f8' }"

    input:
    tuple val(meta),
        path(faa),
        path(blast_raw, stageAs: "in_blast_raw/*"),
        path(blast_best),
        path(blast_unfiltered, stageAs: "in_blast_unfiltered/*"),
        path(blast_filtered, stageAs: "in_blast_filtered/*"),
        path(blast_filtered_best, stageAs: "in_blast_filtered_best/*"),
        path(hmm_compiled, stageAs: "in_hmm_compiled.tsv"),
        path(hmm_per_profile_tsv, stageAs: "in_hmm_tsv/*"),
        path(hmm_per_profile_txt, stageAs: "in_hmm_txt/*"),
        path(hmm_best),
        path(hmm_unfiltered),
        path(hmm_filtered),
        path(merged_sp),
        path(sp_no_fp),
        path(sp_faa),
        path(fp_all, stageAs: "in_fp_compiled.tsv"),
        path(fp_hits),
        path(fp_per_profile_tsv, stageAs: "in_fp_tsv/*"),
        path(fp_per_profile_txt, stageAs: "in_fp_txt/*"),
        path(sp_reannotated),
        path(final_sp),
        path(me_log)
    path mode_file

    output:
    tuple val(meta), path("tmp_intermediate_files.tar.gz"), emit: archive
    tuple val(meta), path("param.conf.gz"), emit: param_conf
    tuple val("${task.process}"), val("tar"), eval('tar --version 2>&1 | head -1 | sed "s/tar (GNU tar) //; s/ Copyright.*//"'), topic: versions, emit: versions_tar
    tuple val("${task.process}"), val("gzip"), eval('gzip --version 2>&1 | head -1 | sed "s/gzip //"'), topic: versions, emit: versions_gzip

    script:
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
    cp in_blast_raw/* ${R}/icescreen_detection_SP/Blast_mode/blastp_output/
    cp in_blast_filtered/* ${R}/icescreen_detection_SP/Blast_mode/filtered_results/
    cp in_blast_filtered_best/* ${R}/icescreen_detection_SP/Blast_mode/filtered_results/
    cp in_blast_unfiltered/* ${R}/icescreen_detection_SP/Blast_mode/unfiltered_results/

    # HMM_mode/
    cp ${hmm_best} ${R}/icescreen_detection_SP/HMM_mode/
    cp in_hmm_compiled.tsv ${R}/icescreen_detection_SP/HMM_mode/hmmscan_output/${S}.tsv
    cp in_hmm_tsv/* ${R}/icescreen_detection_SP/HMM_mode/hmmscan_output/
    cp in_hmm_txt/* ${R}/icescreen_detection_SP/HMM_mode/hmmscan_output/
    cp ${hmm_unfiltered} ${R}/icescreen_detection_SP/HMM_mode/unfiltered_results/
    cp ${hmm_filtered} ${R}/icescreen_detection_SP/HMM_mode/filtered_results/

    # hits_cleaning/
    cp ${merged_sp} ${R}/icescreen_detection_SP/hits_cleaning/
    cp ${sp_faa} ${R}/icescreen_detection_SP/hits_cleaning/
    cp ${sp_no_fp} ${R}/icescreen_detection_SP/hits_cleaning/
    cp ${sp_reannotated} ${R}/icescreen_detection_SP/hits_cleaning/
    cp ${fp_hits} ${R}/icescreen_detection_SP/hits_cleaning/proteins_to_remove/
    cp in_fp_compiled.tsv ${R}/icescreen_detection_SP/hits_cleaning/proteins_to_remove/hmmscan_compiled_output/${S}.tsv
    cp in_fp_tsv/* ${R}/icescreen_detection_SP/hits_cleaning/proteins_to_remove/hmmscan_output/
    cp in_fp_txt/* ${R}/icescreen_detection_SP/hits_cleaning/proteins_to_remove/hmmscan_output/

    tar -czf tmp_intermediate_files.tar.gz results/

    # Generate param.conf: ME detection log params + mode YAML content
    tail -n +8 ${me_log} > param.conf
    tail -n +2 ${mode_file} >> param.conf
    gzip -f param.conf
    """

    stub:
    """
    touch tmp_intermediate_files.tar.gz
    echo "" | gzip > param.conf.gz
    """
}
