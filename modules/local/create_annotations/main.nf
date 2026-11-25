/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CREATE_ANNOTATIONS process
    Generates annotated GFF3, EMBL, and GenBank output files.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process CREATE_ANNOTATIONS {
    tag "${meta.id}"
    label "process_low"

    conda "${moduleDir}/../environment.yml"
    container "ghcr.io/andreassag/icescreen:v1.3.3"

    input:
    tuple val(meta), path(final_sp), path(detected_me), path(genbank)

    output:
    tuple val(meta), path("${prefix}_source.fa.gz"), emit: source_fasta
    tuple val(meta), path("${prefix}_source.gff.gz"), emit: source_gff
    tuple val(meta), path("${prefix}_icescreen.gff.gz"), emit: annotated_gff
    tuple val(meta), path("${prefix}_icescreen.embl.gz"), emit: annotated_embl
    tuple val(meta), path("${prefix}_icescreen.gb.gz"), emit: annotated_gb
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

    # Extract genome sequence from GenBank into FASTA and GFF3
    genbank_to_gff3_and_fasta \\
        -i "\${GB_IN}" \\
        --gff ${prefix}_source.gff \\
        --fasta ${prefix}_source.fa

    # Generate annotated output files (GFF3, EMBL, GenBank)
    generate_annotation_files \\
        ${args} \\
        -s ${final_sp} \\
        -m ${detected_me} \\
        -g "\${GB_IN}" \\
        -o ${prefix}_icescreen

    # Compress outputs
    gzip -f ${prefix}_source.fa
    gzip -f ${prefix}_source.gff
    gzip -f ${prefix}_icescreen.gff
    gzip -f ${prefix}_icescreen.embl
    gzip -f ${prefix}_icescreen.gb
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${prefix}_source.fa.gz
    echo "" | gzip > ${prefix}_source.gff.gz
    echo "" | gzip > ${prefix}_icescreen.gff.gz
    echo "" | gzip > ${prefix}_icescreen.embl.gz
    echo "" | gzip > ${prefix}_icescreen.gb.gz
    """
}
