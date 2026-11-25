/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    QUAST process
    Computes genome assembly contiguity statistics (N50, contig count, total length).
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process QUAST {
    tag "${meta.id}"
    label "process_low"

    conda "bioconda::quast=5.2.0"
    container "quay.io/biocontainers/quast:5.2.0--py39pl5321h4e691d4_3"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("${prefix}/transposed_report.tsv"), emit: transposed_report
    tuple val(meta), path("${prefix}/report.tsv"), emit: report
    tuple val("${task.process}"), val("quast"), eval("quast --version 2>&1 | sed -n '1{s/QUAST v//; p}' || echo 'unknown'"), topic: versions, emit: versions_quast


    script:
    def args = task.ext.args ?: ""
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    quast.py \\
        ${args} \\
        --threads ${task.cpus} \\
        --output-dir ${prefix} \\
        --labels ${meta.id} \\
        ${fasta}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}
    printf 'Assembly\\t# contigs\\tTotal length\\tN50\\n%s\\t1\\t1000000\\t1000000\\n' "${meta.id}" > ${prefix}/transposed_report.tsv
    printf 'Assembly\\t%s\\n# contigs\\t1\\nTotal length\\t1000000\\nN50\\t1000000\\n' "${meta.id}" > ${prefix}/report.tsv
    """
}
