/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ADVENA v2.0.0 Top-Level Workflow
    Orchestrates database provisioning, input normalization & pre-QC, core ICEscreen
    detection, always-on element sequence extraction, opt-in mobile cargo profiling
    (AMRFinderPlus + CARD RGI + DEREP_AMR_CARD, eggNOG-mapper OGs, DefenseFinder,
    antiSMASH-lite), complementary MGE context (MOB-suite, geNomad, ISEScan,
    IntegronFinder, attL/attR), and cohort comparative genomics (skani, MMseqs2,
    MultiQC, pyGenomeViz).
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { samplesheetToList    } from "plugin/nf-schema"

include { PREPARE_DATABASES    } from "../subworkflows/local/prepare_databases/main"
include { INPUT_AND_QC         } from "../subworkflows/local/input_and_qc/main"
include { ICESCREEN            } from "./icescreen"
include { EXTRACT_ME_SEQUENCES } from "../modules/local/extract_me_sequences/main"
include { CARGO_PROFILING      } from "../subworkflows/local/cargo_profiling/main"
include { MGE_CONTEXT          } from "../subworkflows/local/mge_context/main"
include { COHORT_COMPARATIVE   } from "../subworkflows/local/cohort_comparative/main"

def sanitizeVersionValue(raw_version) {
    def text = raw_version == null ? "" : raw_version.toString()
    def first_line = text.readLines()
        .collect { line -> line.replaceAll(/\p{Cntrl}/, " ").trim() }
        .find { line -> line } ?: "unknown"
    return first_line.replace("\\", "\\\\").replace('"', '\\"')
}

def sanitizeYamlKey(raw_key, fallback) {
    def text = raw_key == null ? "" : raw_key.toString()
    def first_line = text.readLines()
        .collect { line -> line.replaceAll(/\p{Cntrl}/, " ").trim() }
        .find { line -> line } ?: fallback
    def safe_key = first_line.replaceAll(/[^A-Za-z0-9._-]/, "_")
    return safe_key ?: fallback
}

def formatVersionEntry(proc_name, tool, raw_version) {
    def simple_name = sanitizeYamlKey(proc_name.toString().tokenize(":").last(), "UNKNOWN_PROCESS")
    def clean_tool = sanitizeYamlKey(tool, "tool")
    def clean_ver = sanitizeVersionValue(raw_version)
    return [simple_name, "${clean_tool}: \"${clean_ver}\""]
}

def formatVersionBlock(simple_name, tool_entries) {
    def body = tool_entries.unique().sort().collect { entry -> "  ${entry}" }.join("\n")
    return "${simple_name}:\n${body}\n"
}

def resolveSamplePath(raw_val, samplesheet_dir, launch_dir, project_dir) {
    if (!raw_val || (raw_val instanceof List && raw_val.isEmpty())) {
        return null
    }
    def p = file(raw_val)
    if (!p.exists()) {
        def abs_p = p.toAbsolutePath()
        def resolved = null
        [launch_dir, project_dir].each { base ->
            if (!resolved && abs_p.startsWith(base)) {
                def candidate = samplesheet_dir.resolve(base.relativize(abs_p).toString())
                if (candidate.exists()) {
                    resolved = candidate
                }
            }
        }
        p = resolved ?: file(raw_val, checkIfExists: true)
    }
    return p
}

workflow ADVENA {

    main:
    if (!params.input) {
        error "ERROR: 'input' parameter is required. Please provide a samplesheet CSV."
    }

    //
    // Resolve opt-in module toggles (--run_all convenience preset)
    //
    def do_qc             = params.run_all || params.run_qc
    def do_taxonomy       = params.run_all || params.run_taxonomy_check
    def do_amr            = params.run_all || params.run_amr
    def do_eggnog         = params.run_all || params.run_eggnog
    def do_defense        = params.run_all || params.run_defense
    def do_bgc            = params.run_all || params.run_bgc
    def do_mge_context    = params.run_all || params.run_mge_context
    def do_att_sites      = params.run_all || params.run_att_sites
    def do_cohort_cluster = params.run_all || params.run_cohort_cluster
    def do_multiqc        = params.run_all || params.run_multiqc
    def do_synteny        = params.run_all || params.run_synteny

    def do_cargo          = do_amr || do_eggnog || do_defense || do_bgc
    def do_cohort_summary = do_qc || do_taxonomy || do_cargo || do_mge_context || do_att_sites || do_cohort_cluster || do_multiqc || do_synteny
    def need_fasta        = do_qc || do_taxonomy || do_mge_context

    //
    // Parse and resolve samplesheet paths (genbank and/or fasta)
    //
    def samplesheet_file = file(params.input, checkIfExists: true)
    def samplesheet_dir = samplesheet_file.toAbsolutePath().parent
    def launch_dir = java.nio.file.Paths.get(workflow.launchDir.toString()).toAbsolutePath()
    def project_dir = java.nio.file.Paths.get(projectDir.toString()).toAbsolutePath()

    def raw_rows = samplesheetToList(params.input, "${projectDir}/assets/schema_input.json")
    def parsed_rows = raw_rows.collect { row ->
        def meta = row[0]
        def raw_gb = row.size() > 1 ? row[1] : null
        def raw_fa = row.size() > 2 ? row[2] : null
        def gb_path = resolveSamplePath(raw_gb, samplesheet_dir, launch_dir, project_dir)
        def fa_path = resolveSamplePath(raw_fa, samplesheet_dir, launch_dir, project_dir)
        [meta, gb_path, fa_path]
    }

    def need_bakta = parsed_rows.any { row -> row[1] == null && row[2] != null }
    ch_samplesheet = channel.fromList(parsed_rows)

    //
    // PHASE 1A: Prepare core ICEscreen + enabled optional databases
    //
    PREPARE_DATABASES(
        need_bakta,
        do_qc,
        do_taxonomy,
        do_amr,
        do_eggnog,
        do_defense,
        do_bgc,
        do_mge_context
    )

    //
    // PHASE 1B: Input normalization (Bakta / FASTA extraction) & optional pre-QC / taxonomy gate
    //
    INPUT_AND_QC(
        ch_samplesheet,
        PREPARE_DATABASES.out.bakta_db,
        PREPARE_DATABASES.out.checkm2_db,
        PREPARE_DATABASES.out.gtdbtk_db,
        do_qc,
        do_taxonomy,
        need_fasta
    )

    //
    // CORE WORKFLOW: ICEscreen detection + always-on element sequence extraction
    //
    ICESCREEN(
        INPUT_AND_QC.out.genbank,
        PREPARE_DATABASES.out.mode_file,
        PREPARE_DATABASES.out.ice_finder_db,
        PREPARE_DATABASES.out.me_conf,
        PREPARE_DATABASES.out.blast_dbs,
        PREPARE_DATABASES.out.sp_hmms,
        PREPARE_DATABASES.out.fp_hmms
    )

    ch_extract_in = ICESCREEN.out.detected_me
        .join(ICESCREEN.out.sp_with_me_ids, by: 0)
        .join(ICESCREEN.out.annotated_gb, by: 0)

    EXTRACT_ME_SEQUENCES(ch_extract_in)

    //
    // PHASE 2: Opt-in mobile cargo profiling (AMR/CARD + DEREP, eggNOG OGs, DefenseFinder, BGC)
    //
    CARGO_PROFILING(
        EXTRACT_ME_SEQUENCES.out.cargo_fna,
        EXTRACT_ME_SEQUENCES.out.bed,
        EXTRACT_ME_SEQUENCES.out.faa,
        EXTRACT_ME_SEQUENCES.out.gff3,
        EXTRACT_ME_SEQUENCES.out.manifest,
        EXTRACT_ME_SEQUENCES.out.me_count,
        ICESCREEN.out.detected_me,
        ICESCREEN.out.annotated_gb,
        ICESCREEN.out.annotated_gff,
        PREPARE_DATABASES.out.amrfinder_db,
        PREPARE_DATABASES.out.card_db,
        PREPARE_DATABASES.out.eggnog_db,
        PREPARE_DATABASES.out.defensefinder_db,
        PREPARE_DATABASES.out.antismash_db,
        do_amr,
        do_eggnog,
        do_defense,
        do_bgc
    )

    //
    // PHASE 3: Opt-in complementary MGE context (MOB-suite, geNomad, ISEScan, IntegronFinder, attL/attR)
    //
    MGE_CONTEXT(
        INPUT_AND_QC.out.fasta,
        ICESCREEN.out.detected_me,
        ICESCREEN.out.annotated_gb,
        ICESCREEN.out.annotated_gff,
        PREPARE_DATABASES.out.mobsuite_db,
        PREPARE_DATABASES.out.genomad_db,
        do_mge_context,
        do_att_sites
    )

    //
    // PHASE 4: Opt-in cohort comparative genomics, MultiQC dashboard, and pyGenomeViz synteny
    //
    ch_versions_topic = channel.topic("versions")
    ch_versions_mqc = channel.empty()
    ch_versions_for_report = ch_versions_topic
    if (do_multiqc) {
        ch_versions_fork = ch_versions_topic.multiMap { item ->
            mqc: item
            all: item
        }
        ch_versions_mqc = ch_versions_fork.mqc
            .distinct()
            .map { proc_name, tool, version -> formatVersionEntry(proc_name, tool, version) }
            .groupTuple(by: 0)
            .map { simple_name, tool_entries -> formatVersionBlock(simple_name, tool_entries) }
            .collectFile(
                name: "software_versions_mqc_versions.yml",
                sort: true,
                newLine: false
            )
        ch_versions_for_report = ch_versions_fork.all
    }

    ch_mqc_extra = INPUT_AND_QC.out.quast_reports
        .mix(INPUT_AND_QC.out.checkm2_reports)
        .mix(INPUT_AND_QC.out.gtdbtk_summaries)
        .mix(INPUT_AND_QC.out.bakta_txts)
        .mix(ch_versions_mqc)

    COHORT_COMPARATIVE(
        ICESCREEN.out.detected_me,
        MGE_CONTEXT.out.contextualized_me,
        CARGO_PROFILING.out.cargo_profile,
        EXTRACT_ME_SEQUENCES.out.bed,
        EXTRACT_ME_SEQUENCES.out.fna,
        EXTRACT_ME_SEQUENCES.out.faa,
        EXTRACT_ME_SEQUENCES.out.manifest,
        EXTRACT_ME_SEQUENCES.out.gbk_slices,
        INPUT_AND_QC.out.qc_summary,
        ch_mqc_extra,
        do_cohort_cluster,
        do_multiqc,
        do_synteny,
        do_cohort_summary
    )

    ch_versions_blocks = ch_versions_for_report
        .mix(COHORT_COMPARATIVE.out.versions_multiqc)
        .distinct()
        .map { proc_name, tool, version -> formatVersionEntry(proc_name, tool, version) }
        .groupTuple(by: 0)
        .map { simple_name, tool_entries -> formatVersionBlock(simple_name, tool_entries) }

    emit:
    detected_me       = ICESCREEN.out.detected_me
    detected_sp       = ICESCREEN.out.detected_sp
    annotations       = ICESCREEN.out.annotated_gff
    me_manifest       = EXTRACT_ME_SEQUENCES.out.manifest
    cargo_profile     = CARGO_PROFILING.out.cargo_profile
    contextualized_me = MGE_CONTEXT.out.contextualized_me
    cohort_summary    = COHORT_COMPARATIVE.out.cohort_summary
    multiqc_report    = COHORT_COMPARATIVE.out.multiqc_report
    versions_multiqc  = COHORT_COMPARATIVE.out.versions_multiqc
    versions_blocks   = ch_versions_blocks
}
