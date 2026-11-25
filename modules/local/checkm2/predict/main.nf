/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CHECKM2_PREDICT process
    Predicts bacterial genome completeness and contamination metrics.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process CHECKM2_PREDICT {
    tag "${meta.id}"
    label "process_medium"

    conda "bioconda::checkm2=1.0.2"
    container "quay.io/biocontainers/checkm2:1.0.2--pyh7cba7a3_0"

    input:
    tuple val(meta), path(fasta)
    path db

    output:
    tuple val(meta), path("${prefix}/quality_report.tsv"), emit: quality_report
    tuple val("${task.process}"), val("checkm2"), eval("checkm2 --version 2>&1 | sed -n '1p' || echo 'unknown'"), topic: versions, emit: versions_checkm2


    script:
    def args = task.ext.args ?: ""
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    export HOME="\$PWD"
    DMND_PATH=\$(find -L ${db} -name "*.dmnd" | sed -n '1p')
    DMND_PATH=\${DMND_PATH:-${db}}

    mkdir -p input_bins
    gzip -cdf ${fasta} > input_bins/${meta.id}.fna

    checkm2 predict \\
        ${args} \\
        --threads ${task.cpus} \\
        --input input_bins \\
        --extension .fna \\
        --database_path "\${DMND_PATH}" \\
        --output-directory ${prefix}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}
    printf 'Name\\tCompleteness\\tContamination\\n%s\\t99.5\\t0.5\\n' "${meta.id}" > ${prefix}/quality_report.tsv
    """
}
