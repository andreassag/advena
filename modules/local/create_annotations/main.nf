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
    tuple val(meta), path("${meta.id}_source.fa.gz"), emit: source_fasta
    tuple val(meta), path("${meta.id}_source.gff.gz"), emit: source_gff
    tuple val(meta), path("${meta.id}_icescreen.gff.gz"), emit: annotated_gff
    tuple val(meta), path("${meta.id}_icescreen.embl.gz"), emit: annotated_embl
    tuple val(meta), path("${meta.id}_icescreen.gb.gz"), emit: annotated_gb
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    # Extract genome sequence from GenBank into FASTA and GFF3
    genbank_to_gff3_and_fasta \\
        -i ${genbank} \\
        --gff ${meta.id}_source.gff \\
        --fasta ${meta.id}_source.fa

    # Generate annotated output files (GFF3, EMBL, GenBank)
    generate_annotation_files \\
        -s ${final_sp} \\
        -m ${detected_me} \\
        -g ${genbank} \\
        -o ${meta.id}_icescreen

    # Compress outputs
    gzip -f ${meta.id}_source.fa
    gzip -f ${meta.id}_source.gff
    gzip -f ${meta.id}_icescreen.gff
    gzip -f ${meta.id}_icescreen.embl
    gzip -f ${meta.id}_icescreen.gb

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version 2>&1 | sed 's/Python //')
        biopython: \$(python3 -c "import Bio; print(Bio.__version__)")
    END_VERSIONS
    """

    stub:
    """
    echo "" | gzip > ${meta.id}_source.fa.gz
    echo "" | gzip > ${meta.id}_source.gff.gz
    echo "" | gzip > ${meta.id}_icescreen.gff.gz
    echo "" | gzip > ${meta.id}_icescreen.embl.gz
    echo "" | gzip > ${meta.id}_icescreen.gb.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: 3.11.0
        biopython: 1.85
    END_VERSIONS
    """
}
