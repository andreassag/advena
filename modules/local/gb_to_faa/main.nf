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
    container "ghcr.io/exterex/icescreen:v1.3.3"

    input:
    tuple val(meta), path(genbank)
    path db

    output:
    tuple val(meta), path("${meta.id}.faa"), emit: faa
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def mode_file = "${db}/mode/${params.phylum}.yml"
    """
    gb_to_faa \\
        -i ${genbank} \\
        -o ${meta.id}.faa \\
        -c ${mode_file}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version 2>&1 | sed 's/Python //')
        biopython: \$(python3 -c "import Bio; print(Bio.__version__)")
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.faa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: 3.11.0
        biopython: 1.85
    END_VERSIONS
    """
}
