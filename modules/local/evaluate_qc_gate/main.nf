/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    EVALUATE_QC_GATE process
    Collates CheckM2, QUAST, and GTDB-Tk outputs into qc_gate_summary.tsv and passed list.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process EVALUATE_QC_GATE {
    tag "qc_gate"
    label "process_single"

    conda "${moduleDir}/../environment.yml"
    container "ghcr.io/andreassag/icescreen:v1.3.3"

    input:
    val sample_ids
    path checkm2_reports, stageAs: "checkm2/?/*"
    path quast_reports, stageAs: "quast/?/*"
    path gtdbtk_summaries, stageAs: "gtdbtk/?/*"

    output:
    path "qc_gate_summary.tsv", emit: summary
    path "passed_samples.txt", emit: passed_samples
    tuple val("${task.process}"), val("python"), eval("python3 --version 2>&1 | sed -n '1{s/Python //; p}'"), topic: versions, emit: versions_python


    script:
    def args = task.ext.args ?: ""
    def samples_arg = sample_ids.collect { s -> "'${s}'" }.join(" ")
    def cm2_arg = checkm2_reports ? "--checkm2 ${checkm2_reports}" : ""
    def quast_arg = quast_reports ? "--quast ${quast_reports}" : ""
    def gtdb_arg = gtdbtk_summaries ? "--gtdbtk ${gtdbtk_summaries}" : ""
    """
    evaluate_qc_gate.py \\
        --samples ${samples_arg} \\
        ${cm2_arg} \\
        ${quast_arg} \\
        ${gtdb_arg} \\
        ${args} \\
        --out-summary qc_gate_summary.tsv \\
        --out-passed passed_samples.txt
    """

    stub:
    def samples_lines = sample_ids.join("\\n")
    """
    printf 'sample\\tcompleteness\\tcontamination\\tn_contigs\\ttotal_length\\tn50\\tgtdb_classification\\tis_bacillota\\tqc_status\\tincluded_downstream\\n' > qc_gate_summary.tsv
    printf '${samples_lines}\\n' > passed_samples.txt
    for s in ${sample_ids.collect { s -> "'${s}'" }.join(" ")}; do
        printf '%s\\t99.50\\t0.50\\t1\\t1000000\\t1000000\\td__Bacteria;p__Bacillota\\tTRUE\\tPASS\\tTRUE\\n' "\$s" >> qc_gate_summary.tsv
    done
    """
}
