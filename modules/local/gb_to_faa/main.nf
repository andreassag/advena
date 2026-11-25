/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    GB_TO_FAA process
    Extracts CDS protein sequences from GenBank files.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process GB_TO_FAA {
    tag "${meta.id}"
    label "process_low"

    conda "${moduleDir}/../environment.yml"
    container "ghcr.io/andreassag/icescreen:v1.3.3"

    input:
    tuple val(meta), path(genbank)
    path mode_file

    output:
    tuple val(meta), path("${prefix}.faa"), emit: faa
    tuple val("${task.process}"), val("icescreen"), eval("python3 -c 'import importlib.metadata; print(importlib.metadata.version(\"icescreen\"))' 2>/dev/null || echo 'unknown'"), topic: versions, emit: versions_icescreen
    tuple val("${task.process}"), val("biopython"), eval("python3 -c 'import Bio; print(Bio.__version__)' 2>/dev/null || echo 'unknown'"), topic: versions, emit: versions_biopython

    script:
    def args = task.ext.args ?: ""
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    GB_IN="${genbank}"
    if [[ "\${GB_IN}" == *.gz ]]; then
        gzip -cdf "${genbank}" > "${prefix}_input.gbff"
        GB_IN="${prefix}_input.gbff"
    fi

    gb_to_faa \\
        ${args} \\
        -i "\${GB_IN}" \\
        -o ${prefix}.faa \\
        -c ${mode_file}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.faa
    """
}
