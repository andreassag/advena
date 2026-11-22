/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PREPARE_DATABASE subworkflow
    Downloads/extracts the ICEscreen database (or reuses a cached db_dir) and builds
    BLAST protein databases and pressed HMMER profile databases.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { UNTAR                         } from "../../../modules/nf-core/untar/main"
include { BLAST_MAKEBLASTDB             } from "../../../modules/nf-core/blast/makeblastdb/main"
include { HMMER_HMMPRESS as HMMPRESS_SP } from "../../../modules/nf-core/hmmer/hmmpress/main"
include { HMMER_HMMPRESS as HMMPRESS_FP } from "../../../modules/nf-core/hmmer/hmmpress/main"

workflow PREPARE_DATABASE {

    take:
    db_url // string: URL or local path to icescreen_database.tar.gz
    db_dir // string|null: optional local cache/database directory
    phylum // string: taxonomic phylum (e.g. "bacillota")

    main:
    def existing_db = null
    if (db_dir) {
        def candidate_nested = file("${db_dir}/icescreen_database")
        def candidate_direct = file("${db_dir}")
        if (file("${candidate_nested}/icescreen.conf").exists()) {
            existing_db = candidate_nested
        } else if (file("${candidate_direct}/icescreen.conf").exists()) {
            existing_db = candidate_direct
        }
    }

    if (existing_db) {
        ch_db_root = channel.value(existing_db)
    } else {
        ch_archive = channel.value([[id: "icescreen_database"], file(db_url, checkIfExists: true)])
        UNTAR(ch_archive)
        ch_db_root = UNTAR.out.untar.map { _meta, dir -> dir }
    }

    ch_mode_file = ch_db_root
        .map { root -> file("${root}/mode/${phylum}.yml", checkIfExists: true) }

    ch_ice_finder_db = ch_db_root
        .map { root -> file("${root}/blastdb/ICE_Finder.db", checkIfExists: true) }

    ch_me_conf = ch_db_root
        .map { root -> file("${root}/icescreen.conf", checkIfExists: true) }

    ch_blast_db_names = channel.of(
        "Relaxase",
        "Couplingprot",
        "Couplingprot2",
        "VirB4",
        "IntTyr",
        "IntSer",
        "DDE"
    )

    ch_blast_faas = ch_db_root
        .combine(ch_blast_db_names)
        .map { root, db_name ->
            [[id: db_name], file("${root}/blastdb/${db_name}.faa", checkIfExists: true)]
        }

    BLAST_MAKEBLASTDB(ch_blast_faas, [])

    ch_sp_hmm_inputs = ch_db_root
        .flatMap { root ->
            def hmm_files = files("${root}/hmmdb/SP_profiles/*.hmm")
            if (!hmm_files) {
                error "No SP HMM profiles found in ${root}/hmmdb/SP_profiles"
            }
            hmm_files.sort { hmm -> hmm.name }.collect { hmm -> [[id: hmm.baseName], hmm] }
        }

    HMMPRESS_SP(ch_sp_hmm_inputs)

    ch_sp_hmms = ch_sp_hmm_inputs
        .join(HMMPRESS_SP.out.compressed_db, by: 0)
        .flatMap { _meta, hmm, h3_files ->
            [hmm] + (h3_files instanceof List ? h3_files : [h3_files])
        }
        .collect()

    ch_fp_hmm_inputs = ch_db_root
        .flatMap { root ->
            def hmm_files = files("${root}/hmmdb/FP_profiles/*.hmm")
            if (!hmm_files) {
                error "No FP HMM profiles found in ${root}/hmmdb/FP_profiles"
            }
            hmm_files.sort { hmm -> hmm.name }.collect { hmm -> [[id: hmm.baseName], hmm] }
        }

    HMMPRESS_FP(ch_fp_hmm_inputs)

    ch_fp_hmms = ch_fp_hmm_inputs
        .join(HMMPRESS_FP.out.compressed_db, by: 0)
        .flatMap { _meta, hmm, h3_files ->
            [hmm] + (h3_files instanceof List ? h3_files : [h3_files])
        }
        .collect()

    emit:
    db_root       = ch_db_root             // value channel: path to extracted database directory
    mode_file     = ch_mode_file           // value channel: path to mode/<phylum>.yml
    ice_finder_db = ch_ice_finder_db       // value channel: path to blastdb/ICE_Finder.db
    me_conf       = ch_me_conf             // value channel: path to icescreen.conf
    blast_dbs     = BLAST_MAKEBLASTDB.out.db // channel: [ val(meta), path(blast_db_dir) ]
    sp_hmms       = ch_sp_hmms             // value channel: list of SP .hmm + .h3* files
    fp_hmms       = ch_fp_hmms             // value channel: list of FP .hmm + .h3* files
}
