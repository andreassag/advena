/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    MULTIQC process
    Aggregates QC reports and advena custom-content panels into an HTML cohort dashboard.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process MULTIQC {
    tag "multiqc"
    label "process_single"

    conda "bioconda::multiqc=1.27.1"
    container "quay.io/biocontainers/multiqc:1.27.1--pyhdfd78af_0"

    input:
    path multiqc_files, stageAs: "?/*"

    output:
    path "multiqc_report.html", emit: report
    path "multiqc_data", emit: data
    tuple val("${task.process}"), val("multiqc"), eval("multiqc --version 2>&1 | sed -n '1{s/multiqc, version //; p}' || echo 'unknown'"), emit: versions_multiqc


    script:
    def args = task.ext.args ?: ""
    """
    export HOME="\$PWD"
    export MPLBACKEND=Agg
    export MPLCONFIGDIR="\$PWD/.mplconfig"
    export XDG_CACHE_HOME="\$PWD/.cache"
    mkdir -p "\${MPLCONFIGDIR}" "\${XDG_CACHE_HOME}"
    multiqc \\
        --force \\
        ${args} \\
        .
    """

    stub:
    """
    mkdir -p multiqc_data
    printf '<!DOCTYPE html><html><head><title>advena MultiQC Stub</title></head><body><h1>MultiQC</h1></body></html>\\n' > multiqc_report.html
    """
}
