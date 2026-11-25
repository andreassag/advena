/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CARGO_PROFILING subworkflow
    Runs AMRFinderPlus (--plus) + CARD RGI + DEREP_AMR_CARD, eggNOG-mapper + OG parsing,
    DefenseFinder, and antiSMASH-lite + BGC intersection, then merges cargo into a
    unified per-element table and enriched advena GFF3.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { AMRFINDERPLUS_RUN       } from "../../../modules/local/amrfinderplus/run/main"
include { RGI_MAIN                } from "../../../modules/local/rgi/main/main"
include { DEREP_AMR_CARD          } from "../../../modules/local/derep_amr_card/main"
include { EGGNOGMAPPER            } from "../../../modules/local/eggnogmapper/main"
include { PARSE_EGGNOG_OGS        } from "../../../modules/local/parse_eggnog_ogs/main"
include { DEFENSEFINDER_RUN       } from "../../../modules/local/defensefinder/run/main"
include { ANTISMASH_ANTISMASHLITE } from "../../../modules/local/antismash/antismashlite/main"
include { INTERSECT_BGC_ME        } from "../../../modules/local/intersect_bgc_me/main"
include { MERGE_CARGO_ANNOTATIONS } from "../../../modules/local/merge_cargo_annotations/main"

workflow CARGO_PROFILING {

    take:
    ch_me_fna           // channel: [ val(meta), path(cargo_query_fna) ]
    ch_me_bed           // channel: [ val(meta), path(me_regions_bed) ]
    ch_me_faa           // channel: [ val(meta), path(faa) ]
    ch_me_gff3          // channel: [ val(meta), path(gff3) ]
    ch_me_manifest      // channel: [ val(meta), path(manifest) ]
    ch_me_count         // channel: [ val(meta), path(count_txt) ]
    ch_detected_me      // channel: [ val(meta), path(detected_me_tsv) ]
    ch_annotated_gb     // channel: [ val(meta), path(annotated_gb_gz) ]
    ch_annotated_gff    // channel: [ val(meta), path(annotated_gff_gz) ]
    ch_amrfinder_db     // value channel
    ch_card_db          // value channel
    ch_eggnog_db        // value channel
    ch_defensefinder_db // value channel
    ch_antismash_db     // value channel
    run_amr             // boolean
    run_eggnog          // boolean
    run_defense         // boolean
    run_bgc             // boolean

    main:
    //
    // Bypass 0-ME samples cleanly when cargo_scope == "elements"
    //
    ch_active_meta = ch_me_count
        .filter { _meta, count_file ->
            if (params.cargo_scope == "genome") {
                return true
            }
            def n = count_file.text.trim()
            return n.isInteger() && n.toInteger() > 0
        }
        .map { meta, _count_file -> [meta.id, meta] }

    ch_active_faa = ch_active_meta
        .join(ch_me_faa.map { meta, f -> [meta.id, f] }, by: 0)
        .map { _id, meta, f -> [meta, f] }

    ch_active_fna = ch_active_meta
        .join(ch_me_fna.map { meta, f -> [meta.id, f] }, by: 0)
        .map { _id, meta, f -> [meta, f] }

    ch_active_bed = ch_active_meta
        .join(ch_me_bed.map { meta, f -> [meta.id, f] }, by: 0)
        .map { _id, meta, f -> [meta, f] }

    ch_active_gff3 = ch_active_meta
        .join(ch_me_gff3.map { meta, f -> [meta.id, f] }, by: 0)
        .map { _id, meta, f -> [meta, f] }

    ch_active_manifest = ch_active_meta
        .join(ch_me_manifest.map { meta, f -> [meta.id, f] }, by: 0)
        .map { _id, meta, f -> [meta, f] }

    //
    // 1. AMRFinderPlus (--plus) + CARD RGI + native DEREP_AMR_CARD
    //
    ch_amr_derep = channel.empty()
    if (run_amr) {
        ch_amrf_in = ch_active_faa
            .map { meta, faa -> [meta.id, meta, faa] }
            .join(ch_active_fna.map { meta, fna -> [meta.id, fna] }, by: 0)
            .join(ch_active_gff3.map { meta, gff -> [meta.id, gff] }, by: 0)
            .map { _id, meta, faa, fna, gff -> [meta, faa, fna, gff] }

        AMRFINDERPLUS_RUN(ch_amrf_in, ch_amrfinder_db)
        RGI_MAIN(ch_active_faa, ch_card_db)

        ch_derep_in = ch_active_manifest
            .map { meta, man -> [meta.id, meta, man] }
            .join(ch_active_bed.map { meta, b -> [meta.id, b] }, by: 0)
            .join(AMRFINDERPLUS_RUN.out.report.map { meta, r -> [meta.id, r] }, by: 0)
            .join(RGI_MAIN.out.tsv.map { meta, t -> [meta.id, t] }, by: 0)
            .map { _id, meta, man, bed, amr, rgi -> [meta, man, bed, amr, rgi] }

        DEREP_AMR_CARD(ch_derep_in)
        ch_amr_derep = DEREP_AMR_CARD.out.derep_tsv
    }

    //
    // 2. eggNOG-mapper (DIAMOND) + hierarchical OG parser (@2 Bacteria + @1239 Bacillota)
    //
    ch_cargo_ogs = channel.empty()
    if (run_eggnog) {
        EGGNOGMAPPER(ch_active_faa, ch_eggnog_db)

        ch_parse_og_in = EGGNOGMAPPER.out.annotations
            .map { meta, ann -> [meta.id, meta, ann] }
            .join(ch_active_manifest.map { meta, man -> [meta.id, man] }, by: 0)
            .map { _id, meta, ann, man -> [meta, ann, man] }

        PARSE_EGGNOG_OGS(ch_parse_og_in)
        ch_cargo_ogs = PARSE_EGGNOG_OGS.out.cargo_ogs
    }

    //
    // 3. DefenseFinder (anti-phage defense & TA systems)
    //
    ch_def_genes = channel.empty()
    ch_def_systems = channel.empty()
    if (run_defense) {
        DEFENSEFINDER_RUN(ch_active_faa, ch_defensefinder_db)
        ch_def_genes = DEFENSEFINDER_RUN.out.genes
        ch_def_systems = DEFENSEFINDER_RUN.out.systems
    }

    //
    // 4. antiSMASH-lite on full annotated GenBank + BGC/ME intersection
    //
    ch_bgc_overlaps = channel.empty()
    if (run_bgc) {
        ch_active_gb = ch_active_meta
            .join(ch_annotated_gb.map { meta, gb -> [meta.id, gb] }, by: 0)
            .map { _id, meta, gb -> [meta, gb] }

        ANTISMASH_ANTISMASHLITE(ch_active_gb, ch_antismash_db)

        ch_bgc_in = ch_active_meta
            .join(ch_detected_me.map { meta, me -> [meta.id, me] }, by: 0)
            .join(ch_me_manifest.map { meta, man -> [meta.id, man] }, by: 0)
            .join(ANTISMASH_ANTISMASHLITE.out.json.map { meta, j -> [meta.id, j] }, by: 0)
            .map { _id, meta, me, man, j -> [meta, me, man, j] }

        INTERSECT_BGC_ME(ch_bgc_in)
        ch_bgc_overlaps = INTERSECT_BGC_ME.out.overlaps
    }

    //
    // 5. Merge all enabled cargo annotations onto ME CDS manifest + enrich GFF3
    //
    ch_cargo_profile = channel.empty()
    ch_enriched_gff = channel.empty()
    if (run_amr || run_eggnog || run_defense || run_bgc) {
        ch_merge_in = ch_me_manifest
            .map { meta, man -> [meta.id, meta, man] }
            .join(ch_annotated_gff.map { meta, g -> [meta.id, g] }, by: 0, remainder: true)
            .join(ch_amr_derep.map { meta, f -> [meta.id, f] }, by: 0, remainder: true)
            .join(ch_cargo_ogs.map { meta, f -> [meta.id, f] }, by: 0, remainder: true)
            .join(ch_def_genes.map { meta, f -> [meta.id, f] }, by: 0, remainder: true)
            .join(ch_def_systems.map { meta, f -> [meta.id, f] }, by: 0, remainder: true)
            .join(ch_bgc_overlaps.map { meta, f -> [meta.id, f] }, by: 0, remainder: true)
            .filter { row -> row.size() >= 3 && row[1] != null }
            .map { row ->
                def meta = row[1]
                def man  = row[2]
                def gff  = (row.size() > 3 && row[3] != null) ? row[3] : []
                def amr  = (row.size() > 4 && row[4] != null) ? row[4] : []
                def og   = (row.size() > 5 && row[5] != null) ? row[5] : []
                def dg   = (row.size() > 6 && row[6] != null) ? row[6] : []
                def ds   = (row.size() > 7 && row[7] != null) ? row[7] : []
                def bgc  = (row.size() > 8 && row[8] != null) ? row[8] : []
                [meta, man, gff, amr, og, dg, ds, bgc]
            }

        MERGE_CARGO_ANNOTATIONS(ch_merge_in)
        ch_cargo_profile = MERGE_CARGO_ANNOTATIONS.out.cargo_profile
        ch_enriched_gff = MERGE_CARGO_ANNOTATIONS.out.enriched_gff
    }

    emit:
    cargo_profile   = ch_cargo_profile
    enriched_gff    = ch_enriched_gff
    amr_derep       = ch_amr_derep
    cargo_ogs       = ch_cargo_ogs
    defense_systems = ch_def_systems
    defense_genes   = ch_def_genes
    bgc_overlaps    = ch_bgc_overlaps
}
