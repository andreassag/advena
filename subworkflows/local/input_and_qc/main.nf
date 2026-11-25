/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    INPUT_AND_QC subworkflow
    Routes GenBank and/or raw FASTA inputs (running Bakta when GenBank is absent),
    extracts FASTA from GenBank when needed, runs QUAST / CheckM2 / GTDB-Tk, and
    applies warn-by-default (or strict) QC/phylum gating.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { BAKTA_BAKTA          } from "../../../modules/local/bakta/bakta/main"
include { EXTRACT_GENOME_FASTA } from "../../../modules/local/extract_genome_fasta/main"
include { QUAST                } from "../../../modules/local/quast/main"
include { CHECKM2_PREDICT      } from "../../../modules/local/checkm2/predict/main"
include { GTDBTK_CLASSIFYWF    } from "../../../modules/local/gtdbtk/classifywf/main"
include { EVALUATE_QC_GATE     } from "../../../modules/local/evaluate_qc_gate/main"

workflow INPUT_AND_QC {

    take:
    ch_samplesheet        // channel: [ val(meta), path_or_null(genbank), path_or_null(fasta) ]
    ch_bakta_db           // value channel: path to Bakta DB or []
    ch_checkm2_db         // value channel: path to CheckM2 DB or []
    ch_gtdbtk_db          // value channel: path to GTDB-Tk DB or []
    run_qc                // boolean
    run_taxonomy          // boolean
    need_fasta_downstream // boolean

    main:
    //
    // STEP 1: Normalize GenBank input (annotate FASTA-only rows with Bakta)
    //
    ch_gb_branched = ch_samplesheet.branch { meta, gb, fa ->
        has_gb: gb != null
            return [meta, gb]
        fasta_only: fa != null
            return [meta, fa]
    }

    BAKTA_BAKTA(ch_gb_branched.fasta_only, ch_bakta_db)

    ch_all_gb = ch_gb_branched.has_gb.mix(BAKTA_BAKTA.out.gbff)
    ch_bakta_txts = BAKTA_BAKTA.out.txt.map { _meta, f -> f }

    //
    // STEP 2: Normalize FASTA input (only if QC, taxonomy, or MGE context requires FASTA)
    //
    ch_all_fa = channel.empty()
    if (need_fasta_downstream) {
        ch_fa_branched = ch_samplesheet.branch { meta, gb, fa ->
            has_fa: fa != null
                return [meta, fa]
            gb_only: gb != null
                return [meta, gb]
        }

        EXTRACT_GENOME_FASTA(ch_fa_branched.gb_only)
        ch_all_fa = ch_fa_branched.has_fa.mix(EXTRACT_GENOME_FASTA.out.fasta)
    }

    //
    // STEP 3: Optional QC (QUAST + CheckM2) and taxonomy verification (GTDB-Tk)
    //
    ch_quast_reports = channel.empty()
    ch_checkm2_reports = channel.empty()
    ch_gtdbtk_summaries = channel.empty()
    ch_qc_summary = channel.empty()

    if (run_qc) {
        QUAST(ch_all_fa)
        CHECKM2_PREDICT(ch_all_fa, ch_checkm2_db)
        ch_quast_reports = QUAST.out.transposed_report.map { _meta, f -> f }
        ch_checkm2_reports = CHECKM2_PREDICT.out.quality_report.map { _meta, f -> f }
    }

    if (run_taxonomy) {
        GTDBTK_CLASSIFYWF(ch_all_fa, ch_gtdbtk_db)
        ch_gtdbtk_summaries = GTDBTK_CLASSIFYWF.out.summary.map { _meta, f -> f }
    }

    ch_filtered_gb = ch_all_gb
    ch_filtered_fa = ch_all_fa

    if (run_qc || run_taxonomy) {
        ch_sample_ids = ch_samplesheet.map { meta, _gb, _fa -> meta.id }.collect()

        EVALUATE_QC_GATE(
            ch_sample_ids,
            ch_checkm2_reports.collect().ifEmpty([]),
            ch_quast_reports.collect().ifEmpty([]),
            ch_gtdbtk_summaries.collect().ifEmpty([])
        )
        ch_qc_summary = EVALUATE_QC_GATE.out.summary

        if (params.strict_qc || params.filter_non_bacillota) {
            ch_passed_set = EVALUATE_QC_GATE.out.passed_samples.map { f ->
                f.readLines().collect { line -> line.trim() }.findAll { line -> line }.toSet()
            }
            ch_filtered_gb = ch_all_gb
                .combine(ch_passed_set)
                .filter { meta, _gb, passed -> passed.contains(meta.id) }
                .map { meta, gb, _passed -> [meta, gb] }

            ch_filtered_fa = ch_all_fa
                .combine(ch_passed_set)
                .filter { meta, _fa, passed -> passed.contains(meta.id) }
                .map { meta, fa, _passed -> [meta, fa] }
        }
    }

    emit:
    genbank          = ch_filtered_gb
    fasta            = ch_filtered_fa
    qc_summary       = ch_qc_summary
    quast_reports    = ch_quast_reports
    checkm2_reports  = ch_checkm2_reports
    gtdbtk_summaries = ch_gtdbtk_summaries
    bakta_txts       = ch_bakta_txts
}
