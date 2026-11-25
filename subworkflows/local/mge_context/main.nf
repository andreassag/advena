/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    MGE_CONTEXT subworkflow
    Runs MOB-suite recon, geNomad, ISEScan, IntegronFinder, and attL/attR detection,
    then contextualizes ICEscreen ME calls with plasmid/prophage/IS/integron/att flags.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { MOBSUITE_RECON   } from "../../../modules/local/mobsuite/recon/main"
include { GENOMAD_ENDTOEND } from "../../../modules/local/genomad/endtoend/main"
include { ISESCAN          } from "../../../modules/local/isescan/main"
include { INTEGRONFINDER   } from "../../../modules/local/integronfinder/main"
include { DETECT_ATT_SITES } from "../../../modules/local/detect_att_sites/main"
include { CONTEXTUALIZE_ME } from "../../../modules/local/contextualize_me/main"

workflow MGE_CONTEXT {

    take:
    ch_fasta         // channel: [ val(meta), path(fasta) ]
    ch_detected_me   // channel: [ val(meta), path(detected_me_tsv) ]
    ch_annotated_gb  // channel: [ val(meta), path(annotated_gb_gz) ]
    ch_annotated_gff // channel: [ val(meta), path(annotated_gff_gz) ]
    ch_mobsuite_db   // value channel
    ch_genomad_db    // value channel
    run_mge_context  // boolean
    run_att_sites    // boolean

    main:
    ch_mob_report = channel.empty()
    ch_genomad_plas = channel.empty()
    ch_genomad_vir = channel.empty()
    ch_isescan = channel.empty()
    ch_integrons = channel.empty()
    ch_att_sites = channel.empty()

    if (run_mge_context) {
        MOBSUITE_RECON(ch_fasta, ch_mobsuite_db)
        GENOMAD_ENDTOEND(ch_fasta, ch_genomad_db)
        ISESCAN(ch_fasta)
        INTEGRONFINDER(ch_fasta)

        ch_mob_report = MOBSUITE_RECON.out.contig_report
        ch_genomad_plas = GENOMAD_ENDTOEND.out.plasmid_summary
        ch_genomad_vir = GENOMAD_ENDTOEND.out.virus_summary
        ch_isescan = ISESCAN.out.results
        ch_integrons = INTEGRONFINDER.out.integrons
    }

    if (run_att_sites) {
        ch_att_in = ch_detected_me
            .map { meta, me -> [meta.id, meta, me] }
            .join(ch_annotated_gb.map { meta, gb -> [meta.id, gb] }, by: 0)
            .join(ch_annotated_gff.map { meta, gff -> [meta.id, gff] }, by: 0, remainder: true)
            .filter { row -> row.size() >= 4 && row[1] != null }
            .map { row ->
                def gff = (row.size() > 4 && row[4] != null) ? row[4] : []
                [row[1], row[2], row[3], gff]
            }

        DETECT_ATT_SITES(ch_att_in)
        ch_att_sites = DETECT_ATT_SITES.out.att_sites
    }

    ch_contextualized_me = channel.empty()
    if (run_mge_context || run_att_sites) {
        ch_ctx_in = ch_detected_me
            .map { meta, me -> [meta.id, meta, me] }
            .join(ch_mob_report.map { meta, f -> [meta.id, f] }, by: 0, remainder: true)
            .join(ch_genomad_plas.map { meta, f -> [meta.id, f] }, by: 0, remainder: true)
            .join(ch_genomad_vir.map { meta, f -> [meta.id, f] }, by: 0, remainder: true)
            .join(ch_isescan.map { meta, f -> [meta.id, f] }, by: 0, remainder: true)
            .join(ch_integrons.map { meta, f -> [meta.id, f] }, by: 0, remainder: true)
            .join(ch_att_sites.map { meta, f -> [meta.id, f] }, by: 0, remainder: true)
            .filter { row -> row.size() >= 3 && row[1] != null }
            .map { row ->
                def meta  = row[1]
                def me    = row[2]
                def m     = (row.size() > 3 && row[3] != null) ? row[3] : []
                def gp    = (row.size() > 4 && row[4] != null) ? row[4] : []
                def gv    = (row.size() > 5 && row[5] != null) ? row[5] : []
                def is_f  = (row.size() > 6 && row[6] != null) ? row[6] : []
                def int_f = (row.size() > 7 && row[7] != null) ? row[7] : []
                def att   = (row.size() > 8 && row[8] != null) ? row[8] : []
                [meta, me, m, gp, gv, is_f, int_f, att]
            }

        CONTEXTUALIZE_ME(ch_ctx_in)
        ch_contextualized_me = CONTEXTUALIZE_ME.out.contextualized_me
    }

    emit:
    contextualized_me = ch_contextualized_me
    att_sites         = ch_att_sites
}
