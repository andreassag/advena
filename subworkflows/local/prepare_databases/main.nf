/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PREPARE_DATABASES subworkflow
    Wraps ICEscreen PREPARE_DATABASE and conditionally prepares databases for enabled
    v2.0.0 optional modules with unified --db_cache_dir and >30 GB safety guards.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { PREPARE_DATABASE         } from "../prepare_database/main"
include { BAKTA_BAKTADBDOWNLOAD    } from "../../../modules/local/bakta/baktadbdownload/main"
include { CHECKM2_DATABASEDOWNLOAD } from "../../../modules/local/checkm2/databasedownload/main"
include { GTDBTK_DOWNLOAD          } from "../../../modules/local/gtdbtk/download/main"
include { AMRFINDERPLUS_UPDATE     } from "../../../modules/local/amrfinderplus/update/main"
include { CARD_DOWNLOAD            } from "../../../modules/local/card_download/main"
include { EGGNOG_DOWNLOAD          } from "../../../modules/local/eggnogmapper/download/main"
include { DEFENSEFINDER_UPDATE     } from "../../../modules/local/defensefinder/update/main"
include { GENOMAD_DOWNLOAD         } from "../../../modules/local/genomad/download/main"
include { MOBSUITE_INIT            } from "../../../modules/local/mobsuite_init/main"

workflow PREPARE_DATABASES {

    take:
    need_bakta      // boolean
    run_qc          // boolean
    run_taxonomy    // boolean
    run_amr         // boolean
    run_eggnog      // boolean
    run_defense     // boolean
    run_bgc         // boolean
    run_mge_context // boolean

    main:
    //
    // Core ICEscreen database (always enabled)
    //
    def icescreen_db_dir = params.db_dir ?: (params.db_cache_dir ? "${params.db_cache_dir}/icescreen" : null)
    PREPARE_DATABASE(
        params.db_url,
        icescreen_db_dir,
        params.phylum
    )

    //
    // Large DB (>30 GB) startup safety guards when not in -stub-run
    //
    if (!workflow.stubRun) {
        if (run_taxonomy && !params.gtdbtk_db && !params.auto_download_large_dbs) {
            error(
                "ERROR: --run_taxonomy_check requires a GTDB-Tk reference database (>80 GB).\n" +
                "Provide --gtdbtk_db /path/to/gtdbtk_db (download via https://ecogenomics.github.io/GTDBTk/installing/index.html) " +
                "or explicitly set --auto_download_large_dbs true."
            )
        }
        if (run_eggnog && !params.eggnog_db && !params.auto_download_large_dbs) {
            error(
                "ERROR: --run_eggnog requires an eggNOG-mapper database (~50 GB).\n" +
                "Download via 'download_eggnog_data.py -y --data_dir /path/to/eggnog_db' and pass --eggnog_db /path/to/eggnog_db, " +
                "or explicitly set --auto_download_large_dbs true."
            )
        }
        if (need_bakta && params.bakta_db_type == "full" && !params.bakta_db && !params.auto_download_large_dbs) {
            error(
                "ERROR: --bakta_db_type full requires >35 GB download.\n" +
                "Provide --bakta_db /path/to/bakta_db, set --auto_download_large_dbs true, or use default --bakta_db_type light."
            )
        }
    }

    //
    // Bakta database (only if raw FASTA inputs require annotation)
    //
    ch_bakta_db = channel.value([])
    if (need_bakta) {
        if (params.bakta_db) {
            ch_bakta_db = channel.value(file(params.bakta_db, checkIfExists: true))
        } else {
            BAKTA_BAKTADBDOWNLOAD(params.bakta_db_type)
            ch_bakta_db = BAKTA_BAKTADBDOWNLOAD.out.db
        }
    }

    //
    // CheckM2 database (--run_qc)
    //
    ch_checkm2_db = channel.value([])
    if (run_qc) {
        if (params.checkm2_db) {
            ch_checkm2_db = channel.value(file(params.checkm2_db, checkIfExists: true))
        } else {
            CHECKM2_DATABASEDOWNLOAD()
            ch_checkm2_db = CHECKM2_DATABASEDOWNLOAD.out.db
        }
    }

    //
    // GTDB-Tk database (--run_taxonomy_check)
    //
    ch_gtdbtk_db = channel.value([])
    if (run_taxonomy) {
        if (params.gtdbtk_db) {
            ch_gtdbtk_db = channel.value(file(params.gtdbtk_db, checkIfExists: true))
        } else if (params.auto_download_large_dbs || workflow.stubRun) {
            GTDBTK_DOWNLOAD()
            ch_gtdbtk_db = GTDBTK_DOWNLOAD.out.db
        }
    }

    //
    // AMRFinderPlus & CARD databases (--run_amr)
    //
    ch_amrfinder_db = channel.value([])
    ch_card_db = channel.value([])
    if (run_amr) {
        if (params.amrfinderplus_db) {
            ch_amrfinder_db = channel.value(file(params.amrfinderplus_db, checkIfExists: true))
        } else {
            AMRFINDERPLUS_UPDATE()
            ch_amrfinder_db = AMRFINDERPLUS_UPDATE.out.db
        }

        if (params.card_db) {
            ch_card_db = channel.value(file(params.card_db, checkIfExists: true))
        } else {
            CARD_DOWNLOAD()
            ch_card_db = CARD_DOWNLOAD.out.db
        }
    }

    //
    // eggNOG-mapper database (--run_eggnog)
    //
    ch_eggnog_db = channel.value([])
    if (run_eggnog) {
        if (params.eggnog_db) {
            ch_eggnog_db = channel.value(file(params.eggnog_db, checkIfExists: true))
        } else if (params.auto_download_large_dbs || workflow.stubRun) {
            EGGNOG_DOWNLOAD()
            ch_eggnog_db = EGGNOG_DOWNLOAD.out.db
        }
    }

    //
    // DefenseFinder models (--run_defense)
    //
    ch_defensefinder_db = channel.value([])
    if (run_defense) {
        if (params.defensefinder_db) {
            ch_defensefinder_db = channel.value(file(params.defensefinder_db, checkIfExists: true))
        } else {
            DEFENSEFINDER_UPDATE()
            ch_defensefinder_db = DEFENSEFINDER_UPDATE.out.db
        }
    }

    //
    // antiSMASH database (--run_bgc, optional override)
    //
    ch_antismash_db = channel.value([])
    if (run_bgc && params.antismash_db) {
        ch_antismash_db = channel.value(file(params.antismash_db, checkIfExists: true))
    }

    //
    // MOB-suite & geNomad databases (--run_mge_context)
    //
    ch_mobsuite_db = channel.value([])
    ch_genomad_db = channel.value([])
    if (run_mge_context) {
        if (params.mobsuite_db) {
            ch_mobsuite_db = channel.value(file(params.mobsuite_db, checkIfExists: true))
        } else {
            MOBSUITE_INIT()
            ch_mobsuite_db = MOBSUITE_INIT.out.db
        }

        if (params.genomad_db) {
            ch_genomad_db = channel.value(file(params.genomad_db, checkIfExists: true))
        } else {
            GENOMAD_DOWNLOAD()
            ch_genomad_db = GENOMAD_DOWNLOAD.out.db
        }
    }

    emit:
    db_root          = PREPARE_DATABASE.out.db_root
    mode_file        = PREPARE_DATABASE.out.mode_file
    ice_finder_db    = PREPARE_DATABASE.out.ice_finder_db
    me_conf          = PREPARE_DATABASE.out.me_conf
    blast_dbs        = PREPARE_DATABASE.out.blast_dbs
    sp_hmms          = PREPARE_DATABASE.out.sp_hmms
    fp_hmms          = PREPARE_DATABASE.out.fp_hmms
    bakta_db         = ch_bakta_db
    checkm2_db       = ch_checkm2_db
    gtdbtk_db        = ch_gtdbtk_db
    amrfinder_db     = ch_amrfinder_db
    card_db          = ch_card_db
    eggnog_db        = ch_eggnog_db
    defensefinder_db = ch_defensefinder_db
    antismash_db     = ch_antismash_db
    mobsuite_db      = ch_mobsuite_db
    genomad_db       = ch_genomad_db
}
