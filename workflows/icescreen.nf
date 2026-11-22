/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ICESCREEN workflow
    Detection and annotation of ICEs and IMEs in Bacillota genomes.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { samplesheetToList                 } from "plugin/nf-schema"

include { PREPARE_DATABASE                  } from "../subworkflows/local/prepare_database/main"
include { BLAST_BLASTP                      } from "../modules/nf-core/blast/blastp/main"
include { GB_TO_FAA                         } from "../modules/local/gb_to_faa/main"
include { PROCESS_BLASTP                    } from "../modules/local/process_blastp/main"
include { HMMER_HMMSCAN as HMMSCAN_SP       } from "../modules/local/hmmer/hmmscan/main"
include { PROCESS_HMMSCAN                   } from "../modules/local/process_hmmscan/main"
include { MERGE_SP                          } from "../modules/local/merge_sp/main"
include { EXTRACT_SP_FAA                    } from "../modules/local/extract_sp_faa/main"
include { HMMER_HMMSCAN as HMMSCAN_FP       } from "../modules/local/hmmer/hmmscan/main"
include { PROCESS_HMMSCAN_FP                } from "../modules/local/process_hmmscan_fp/main"
include { REANNOT_SP                        } from "../modules/local/reannot_sp/main"
include { DETECT_ME                         } from "../modules/local/detect_me/main"
include { CREATE_ANNOTATIONS                } from "../modules/local/create_annotations/main"
include { PACKAGE_RESULTS                   } from "../modules/local/package_results/main"

workflow ICESCREEN {

    main:
    //
    // Validate and parse samplesheet via nf-schema
    //
    if (!params.input) {
        error "ERROR: 'input' parameter is required. Please provide a samplesheet CSV."
    }

    def samplesheet_file = file(params.input, checkIfExists: true)
    def samplesheet_dir = samplesheet_file.toAbsolutePath().parent
    def launch_dir = java.nio.file.Paths.get(workflow.launchDir.toString()).toAbsolutePath()
    def project_dir = java.nio.file.Paths.get(projectDir.toString()).toAbsolutePath()

    ch_input = channel
        .fromList(samplesheetToList(params.input, "${projectDir}/assets/schema_input.json"))
        .map { meta, genbank ->
            def gb_path = file(genbank)
            if (!gb_path.exists()) {
                def abs_gb = gb_path.toAbsolutePath()
                def resolved = null
                [launch_dir, project_dir].each { base ->
                    if (!resolved && abs_gb.startsWith(base)) {
                        def candidate = samplesheet_dir.resolve(base.relativize(abs_gb).toString())
                        if (candidate.exists()) {
                            resolved = candidate
                        }
                    }
                }
                gb_path = resolved ?: file(genbank, checkIfExists: true)
            }
            [meta, gb_path]
        }

    //
    // STEP 1: Prepare databases (extract/cache, index BLASTP DBs and press HMMs)
    //
    PREPARE_DATABASE(
        params.db_url,
        params.db_dir,
        params.phylum
    )

    //
    // STEP 2: Extract CDS protein sequences from GenBank
    //
    GB_TO_FAA(ch_input, PREPARE_DATABASE.out.mode_file)

    //
    // STEP 3: BLASTP signature protein search (scattered across 7 databases)
    //
    ch_blast_in = GB_TO_FAA.out.faa
        .combine(PREPARE_DATABASE.out.blast_dbs)
        .multiMap { meta, faa, meta_db, db_dir ->
            query: [meta, faa]
            db: [meta_db, db_dir]
        }

    BLAST_BLASTP(ch_blast_in.query, ch_blast_in.db, "tsv")

    ch_blast_collected = BLAST_BLASTP.out.tsv
        .groupTuple(by: 0, size: 7)

    PROCESS_BLASTP(
        ch_blast_collected,
        PREPARE_DATABASE.out.mode_file,
        PREPARE_DATABASE.out.ice_finder_db
    )

    //
    // STEP 4: HMM signature protein scan and filtering
    //
    HMMSCAN_SP(GB_TO_FAA.out.faa, PREPARE_DATABASE.out.sp_hmms)

    PROCESS_HMMSCAN(HMMSCAN_SP.out.results, PREPARE_DATABASE.out.mode_file)

    //
    // STEP 5: Merge BLAST and HMM best hits
    //
    ch_merge = PROCESS_BLASTP.out.best_hits
        .join(PROCESS_HMMSCAN.out.best_hits, by: 0)

    MERGE_SP(ch_merge)

    //
    // STEP 6: False positive screening (extract candidate SPs, hmmscan FP profiles, filter)
    //
    ch_extract_sp = MERGE_SP.out.merged
        .join(GB_TO_FAA.out.faa, by: 0)

    EXTRACT_SP_FAA(ch_extract_sp)

    HMMSCAN_FP(EXTRACT_SP_FAA.out.sp_faa, PREPARE_DATABASE.out.fp_hmms)

    ch_process_fp = MERGE_SP.out.merged
        .join(HMMSCAN_FP.out.results, by: 0)

    PROCESS_HMMSCAN_FP(ch_process_fp, PREPARE_DATABASE.out.mode_file)

    //
    // STEP 7: Re-annotate signature proteins (XerS classification)
    //
    ch_reannot = PROCESS_HMMSCAN_FP.out.sp_clean
        .join(PROCESS_BLASTP.out.intyr_hits, by: 0)

    REANNOT_SP(ch_reannot, PREPARE_DATABASE.out.mode_file)

    //
    // STEP 8: Detect ICE/IME mobile element structures
    //
    ch_detect = REANNOT_SP.out.final_sp
        .join(ch_input, by: 0)

    DETECT_ME(
        ch_detect,
        PREPARE_DATABASE.out.mode_file,
        PREPARE_DATABASE.out.me_conf
    )

    //
    // STEP 9: Generate annotated output files (GFF3, EMBL, GenBank)
    //
    ch_annotate = REANNOT_SP.out.final_sp
        .join(DETECT_ME.out.detected_me, by: 0)
        .join(ch_input, by: 0)

    CREATE_ANNOTATIONS(ch_annotate)

    //
    // STEP 10: Archive intermediate files and generate run summary (param.conf.gz)
    //
    ch_package = GB_TO_FAA.out.faa
        .join(PROCESS_BLASTP.out.blast_raw, by: 0)
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
        .join(PROCESS_HMMSCAN_FP.out.sp_clean, by: 0)
        .join(EXTRACT_SP_FAA.out.sp_faa, by: 0)
        .join(HMMSCAN_FP.out.results, by: 0)
        .join(PROCESS_HMMSCAN_FP.out.fp_hits, by: 0)
        .join(HMMSCAN_FP.out.per_profile_tsv, by: 0)
        .join(HMMSCAN_FP.out.per_profile_txt, by: 0)
        .join(REANNOT_SP.out.reannotated, by: 0)
        .join(REANNOT_SP.out.final_sp, by: 0)
        .join(DETECT_ME.out.log, by: 0)

    PACKAGE_RESULTS(ch_package, PREPARE_DATABASE.out.mode_file)

    emit:
    detected_me   = DETECT_ME.out.detected_me
    detected_sp   = REANNOT_SP.out.final_sp
    annotations   = CREATE_ANNOTATIONS.out.annotated_gff
}
