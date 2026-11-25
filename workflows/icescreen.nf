/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ICESCREEN workflow
    Core detection and annotation of ICEs and IMEs in Bacillota genomes.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { BLAST_BLASTP                } from "../modules/nf-core/blast/blastp/main"
include { GB_TO_FAA                   } from "../modules/local/gb_to_faa/main"
include { PROCESS_BLASTP              } from "../modules/local/process_blastp/main"
include { HMMER_HMMSCAN as HMMSCAN_SP } from "../modules/local/hmmer/hmmscan/main"
include { PROCESS_HMMSCAN             } from "../modules/local/process_hmmscan/main"
include { MERGE_SP                    } from "../modules/local/merge_sp/main"
include { EXTRACT_SP_FAA              } from "../modules/local/extract_sp_faa/main"
include { HMMER_HMMSCAN as HMMSCAN_FP } from "../modules/local/hmmer/hmmscan/main"
include { PROCESS_HMMSCAN_FP          } from "../modules/local/process_hmmscan_fp/main"
include { REANNOT_SP                  } from "../modules/local/reannot_sp/main"
include { DETECT_ME                   } from "../modules/local/detect_me/main"
include { CREATE_ANNOTATIONS          } from "../modules/local/create_annotations/main"
include { PACKAGE_RESULTS             } from "../modules/local/package_results/main"

workflow ICESCREEN {

    take:
    ch_input      // channel: [ val(meta), path(genbank) ]
    ch_mode_file  // value channel: path(mode_file)
    ch_ice_finder // value channel: path(ice_finder_db)
    ch_me_conf    // value channel: path(me_conf)
    ch_blast_dbs  // channel: [ val(meta_db), path(blast_db_dir) ]
    ch_sp_hmms    // value channel: list of SP .hmm + .h3* files
    ch_fp_hmms    // value channel: list of FP .hmm + .h3* files

    main:
    //
    // STEP 1: Extract CDS protein sequences from GenBank
    //
    GB_TO_FAA(ch_input, ch_mode_file)

    //
    // STEP 2: BLASTP signature protein search (scattered across 7 databases)
    //
    ch_blast_in = GB_TO_FAA.out.faa
        .combine(ch_blast_dbs)
        .multiMap { meta, faa, meta_db, db_dir ->
            query: [meta, faa]
            db: [meta_db, db_dir]
        }

    BLAST_BLASTP(ch_blast_in.query, ch_blast_in.db, "tsv")

    ch_blast_collected = BLAST_BLASTP.out.tsv
        .groupTuple(by: 0, size: 7)

    PROCESS_BLASTP(
        ch_blast_collected,
        ch_mode_file,
        ch_ice_finder
    )

    //
    // STEP 3: HMM signature protein scan and filtering
    //
    HMMSCAN_SP(GB_TO_FAA.out.faa, ch_sp_hmms)

    PROCESS_HMMSCAN(HMMSCAN_SP.out.results, ch_mode_file)

    //
    // STEP 4: Merge BLAST and HMM best hits
    //
    ch_merge = PROCESS_BLASTP.out.best_hits
        .join(PROCESS_HMMSCAN.out.best_hits, by: 0)

    MERGE_SP(ch_merge)

    //
    // STEP 5: False positive screening (extract candidate SPs, hmmscan FP profiles, filter)
    //
    ch_extract_sp = MERGE_SP.out.merged
        .join(GB_TO_FAA.out.faa, by: 0)

    EXTRACT_SP_FAA(ch_extract_sp)

    HMMSCAN_FP(EXTRACT_SP_FAA.out.sp_faa, ch_fp_hmms)

    ch_process_fp = MERGE_SP.out.merged
        .join(HMMSCAN_FP.out.results, by: 0)

    PROCESS_HMMSCAN_FP(ch_process_fp, ch_mode_file)

    //
    // STEP 6: Re-annotate signature proteins (XerS classification)
    //
    ch_reannot = PROCESS_HMMSCAN_FP.out.sp_clean
        .join(PROCESS_BLASTP.out.intyr_hits, by: 0)

    REANNOT_SP(ch_reannot, ch_mode_file)

    //
    // STEP 7: Detect ICE/IME mobile element structures
    //
    ch_detect = REANNOT_SP.out.final_sp
        .join(ch_input, by: 0)

    DETECT_ME(
        ch_detect,
        ch_mode_file,
        ch_me_conf
    )

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

    PACKAGE_RESULTS(ch_package, ch_mode_file)

    emit:
    detected_me    = DETECT_ME.out.detected_me
    sp_with_me_ids = DETECT_ME.out.sp_with_me_ids
    summary        = DETECT_ME.out.summary
    detected_sp    = REANNOT_SP.out.final_sp
    annotated_gff  = CREATE_ANNOTATIONS.out.annotated_gff
    annotated_gb   = CREATE_ANNOTATIONS.out.annotated_gb
    source_fasta   = CREATE_ANNOTATIONS.out.source_fasta
}
