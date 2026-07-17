/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ICESCREEN workflow
    Detection and annotation of ICEs and IMEs in Bacillota genomes.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { DOWNLOAD_DB } from "../modules/local/download_db/main"
include { GB_TO_FAA } from "../modules/local/gb_to_faa/main"
include { BLASTP_SEARCH } from "../modules/local/blastp_search/main"
include { PROCESS_BLASTP } from "../modules/local/process_blastp/main"
include { HMMSCAN_SP } from "../modules/local/hmmscan_sp/main"
include { PROCESS_HMMSCAN } from "../modules/local/process_hmmscan/main"
include { MERGE_SP } from "../modules/local/merge_sp/main"
include { FP_SCREENING } from "../modules/local/fp_screening/main"
include { REANNOT_SP } from "../modules/local/reannot_sp/main"
include { DETECT_ME } from "../modules/local/detect_me/main"
include { CREATE_ANNOTATIONS } from "../modules/local/create_annotations/main"
include { PACKAGE_RESULTS } from "../modules/local/package_results/main"

workflow ICESCREEN {

    //
    // Validate and parse samplesheet
    //
    if (!params.input) {
        error "ERROR: 'input' parameter is required. Please provide a samplesheet CSV."
    }

    channel
        .fromPath(params.input, checkIfExists: true)
        .splitCsv(header: true, strip: true)
        .map { row ->
            def meta = [id: row.sample]
            def samplesheet_dir = file(params.input).parent
            def genbank = row.genbank.startsWith("/")
                ? file(row.genbank, checkIfExists: true)
                : file("${samplesheet_dir}/${row.genbank}", checkIfExists: true)
            [meta, genbank]
        }
        .set { ch_input }

    // Check if a directory was provided and if the database already exists there
    def db_folder = params.db_dir ? file("${params.db_dir}/icescreen_database") : null

    if ( db_folder && db_folder.exists() && db_folder.isDirectory() ) {
        ch_db = channel.value(db_folder)
    } else {
        ch_db = DOWNLOAD_DB()
        ch_db = DOWNLOAD_DB.out.db.first()
    }

    //
    // STEP 1b: Extract CDS protein sequences from GenBank
    //
    GB_TO_FAA(ch_input, ch_db)

    //
    // STEP 2a: BLASTP signature protein search (scattered across 7 databases)
    //
    ch_blast_databases = channel.of(
        "Relaxase",
        "Couplingprot",
        "Couplingprot2",
        "VirB4",
        "IntTyr",
        "IntSer",
        "DDE"
    )

    ch_blast_input = GB_TO_FAA.out.faa
        .combine(ch_blast_databases)

    BLASTP_SEARCH(ch_blast_input, ch_db)

    //
    // STEP 3b: Collect all BLAST results per sample and process
    //
    ch_blast_collected = BLASTP_SEARCH.out.results.groupTuple(by: 0)

    PROCESS_BLASTP(ch_blast_collected, ch_db)

    //
    // STEP 2b: HMM signature protein scan (all profiles in one process)
    //
    HMMSCAN_SP(GB_TO_FAA.out.faa, ch_db)

    //
    // STEP 3b: Filter HMM results
    //
    PROCESS_HMMSCAN(HMMSCAN_SP.out.results, ch_db)

    //
    // STEP 3: Merge BLAST and HMM best hits
    //
    ch_merge = PROCESS_BLASTP.out.best_hits
        .join(PROCESS_HMMSCAN.out.best_hits, by: 0)

    MERGE_SP(ch_merge)

    //
    // STEP 5: False positive screening
    //
    ch_fp = MERGE_SP.out.merged
        .join(GB_TO_FAA.out.faa, by: 0)

    FP_SCREENING(ch_fp, ch_db)

    //
    // STEP 6: Re-annotate signature proteins (XerS classification)
    //
    ch_reannot = FP_SCREENING.out.sp_clean
        .join(PROCESS_BLASTP.out.intyr_hits, by: 0)

    REANNOT_SP(ch_reannot, ch_db)

    //
    // STEP 7: Detect ICE/IME mobile element structures
    //
    ch_detect = REANNOT_SP.out.final_sp
        .join(ch_input, by: 0)

    DETECT_ME(ch_detect, ch_db)

    //
    // STEP 8: Generate annotated output files (GFF3, EMBL, GenBank)
    //
    ch_annotate = REANNOT_SP.out.final_sp
        .join(DETECT_ME.out.detected_me, by: 0)
        .join(ch_input, by: 0)

    CREATE_ANNOTATIONS(ch_annotate)

    //
    // STEP 9: Archive intermediate files and generate run summary (param.conf.gz)
    //
    ch_blast_raw_for_archive = BLASTP_SEARCH.out.results
        .map { meta, db_name, file -> [meta, file] }
        .groupTuple(by: 0)

    ch_package = GB_TO_FAA.out.faa
        .join(ch_blast_raw_for_archive, by: 0)
        .join(PROCESS_BLASTP.out.best_hits, by: 0)
        .join(PROCESS_BLASTP.out.blast_unfiltered, by: 0)
        .join(PROCESS_BLASTP.out.blast_filtered, by: 0)
        .join(PROCESS_BLASTP.out.blast_filtered_best, by: 0)
        .join(HMMSCAN_SP.out.results, by: 0)
        .join(HMMSCAN_SP.out.per_profile_tsv, by: 0)
        .join(HMMSCAN_SP.out.per_profile_txt, by: 0)
        .join(PROCESS_HMMSCAN.out.best_hits, by: 0)
        .join(PROCESS_HMMSCAN.out.hmm_unfiltered, by: 0)
        .join(PROCESS_HMMSCAN.out.hmm_filtered, by: 0)
        .join(MERGE_SP.out.merged, by: 0)
        .join(FP_SCREENING.out.sp_clean, by: 0)
        .join(FP_SCREENING.out.sp_faa, by: 0)
        .join(FP_SCREENING.out.fp_all, by: 0)
        .join(FP_SCREENING.out.fp_hits, by: 0)
        .join(FP_SCREENING.out.fp_per_profile_tsv, by: 0)
        .join(FP_SCREENING.out.fp_per_profile_txt, by: 0)
        .join(REANNOT_SP.out.reannotated, by: 0)
        .join(REANNOT_SP.out.final_sp, by: 0)
        .join(DETECT_ME.out.log, by: 0)

    PACKAGE_RESULTS(ch_package, ch_db)
}
