/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    MERGE_SP process
    Merges BLAST and HMM signature protein results.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process MERGE_SP {
    tag "${meta.id}"
    label "process_low"

    conda "${moduleDir}/../environment.yml"
    container "ghcr.io/andreassag/icescreen:v1.3.3"

    input:
    tuple val(meta), path(blast_best), path(hmm_best)

    output:
    tuple val(meta), path("${prefix}_detected_SP_source.tsv"), emit: merged
    tuple val("${task.process}"), val("icescreen"), eval("python3 -c 'import importlib.metadata; print(importlib.metadata.version(\"icescreen\"))' 2>/dev/null || echo 'unknown'"), topic: versions, emit: versions_icescreen
    tuple val("${task.process}"), val("pandas"), eval("python3 -c 'import pandas; print(pandas.__version__)' 2>/dev/null || echo 'unknown'"), topic: versions, emit: versions_pandas

    script:
    def args = task.ext.args ?: ""
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    merge_SP_results \\
        ${args} \\
        --blastres ${blast_best} \\
        --hmmres ${hmm_best} \\
        -o ${prefix}_detected_SP_source.tsv
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_detected_SP_source.tsv
    """
}
