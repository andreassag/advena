/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COHORT_COMPARATIVE subworkflow
    Cross-genome ICE/IME clustering (skani + MMseqs2), pan-ME eggNOG OG matrix & cohort
    summary, MultiQC dashboard, and pyGenomeViz cluster synteny visualizations.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { SKANI_TRIANGLE      } from "../../../modules/local/skani/triangle/main"
include { MMSEQS_EASYCLUSTER  } from "../../../modules/local/mmseqs/easycluster/main"
include { CLUSTER_MES         } from "../../../modules/local/cluster_mes/main"
include { COHORT_SUMMARY      } from "../../../modules/local/cohort_summary/main"
include { MULTIQC             } from "../../../modules/local/multiqc/main"
include { PYGENOMEVIZ_SYNTENY } from "../../../modules/local/pygenomeviz_synteny/main"

workflow COHORT_COMPARATIVE {

    take:
    ch_detected_me       // channel: [ val(meta), path(detected_me_tsv) ]
    ch_contextualized_me // channel: [ val(meta), path(contextualized_me_tsv) ]
    ch_cargo_profile     // channel: [ val(meta), path(cargo_profile_tsv) ]
    ch_me_bed            // channel: [ val(meta), path(bed) ]
    ch_me_fna            // channel: [ val(meta), path(fna_gz) ]
    ch_me_faa            // channel: [ val(meta), path(faa_gz) ]
    ch_me_manifest       // channel: [ val(meta), path(manifest_tsv) ]
    ch_gbk_slices        // channel: [ val(meta), path(gbk_gz) ]
    ch_qc_summary        // channel: path(qc_gate_summary.tsv)
    ch_mqc_extra         // channel: path(qc / bakta / gtdbtk files)
    run_cohort_cluster   // boolean
    run_multiqc          // boolean
    run_synteny          // boolean
    run_cohort_summary   // boolean

    main:
    ch_clusters_tsv = channel.empty()
    if (run_cohort_cluster || run_synteny) {
        SKANI_TRIANGLE(ch_me_fna.map { _meta, f -> f }.collect().ifEmpty([]))
        MMSEQS_EASYCLUSTER(ch_me_faa.map { _meta, f -> f }.collect().ifEmpty([]))

        CLUSTER_MES(
            ch_me_bed.map { _meta, f -> f }.collect().ifEmpty([]),
            ch_me_manifest.map { _meta, f -> f }.collect().ifEmpty([]),
            SKANI_TRIANGLE.out.matrix,
            MMSEQS_EASYCLUSTER.out.clusters
        )
        ch_clusters_tsv = CLUSTER_MES.out.clusters_tsv
    }

    ch_cohort_summary = channel.empty()
    ch_og_matrix = channel.empty()
    ch_mqc_report = channel.empty()
    ch_versions_multiqc = channel.empty()
    ch_synteny_svgs = channel.empty()

    if (run_cohort_summary || run_cohort_cluster || run_multiqc || run_synteny) {
        COHORT_SUMMARY(
            ch_detected_me.map { _meta, f -> f }.collect().ifEmpty([]),
            ch_contextualized_me.map { _meta, f -> f }.collect().ifEmpty([]),
            ch_cargo_profile.map { _meta, f -> f }.collect().ifEmpty([]),
            ch_qc_summary.ifEmpty([]),
            ch_clusters_tsv.ifEmpty([])
        )
        ch_cohort_summary = COHORT_SUMMARY.out.cohort_summary
        ch_og_matrix = COHORT_SUMMARY.out.og_matrix

        if (run_multiqc) {
            ch_mqc_all = ch_mqc_extra
                .mix(COHORT_SUMMARY.out.mqc_icescreen)
                .mix(COHORT_SUMMARY.out.mqc_amr)
                .mix(COHORT_SUMMARY.out.mqc_cog)
                .mix(COHORT_SUMMARY.out.mqc_mge)
                .collect()

            MULTIQC(ch_mqc_all)
            ch_mqc_report = MULTIQC.out.report
            ch_versions_multiqc = MULTIQC.out.versions_multiqc
        }
    }

    if (run_synteny) {
        ch_all_gbks = ch_gbk_slices
            .flatMap { _meta, gbk_files ->
                if (gbk_files instanceof List) {
                    return gbk_files
                }
                return [gbk_files]
            }
            .collect()
            .ifEmpty([])

        PYGENOMEVIZ_SYNTENY(
            ch_clusters_tsv,
            ch_all_gbks,
            ch_cargo_profile.map { _meta, f -> f }.collect().ifEmpty([])
        )
        ch_synteny_svgs = PYGENOMEVIZ_SYNTENY.out.svgs
    }

    emit:
    clusters_tsv     = ch_clusters_tsv
    cohort_summary   = ch_cohort_summary
    og_matrix        = ch_og_matrix
    multiqc_report   = ch_mqc_report
    versions_multiqc = ch_versions_multiqc
    synteny_svgs     = ch_synteny_svgs
}
