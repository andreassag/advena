/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    EXTRACT_GENOME_FASTA process
    Extracts nucleotide FASTA sequences from GenBank input files for QC and MGE context.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process EXTRACT_GENOME_FASTA {
    tag "${meta.id}"
    label "process_single"

    conda "${moduleDir}/../environment.yml"
    container "ghcr.io/andreassag/icescreen:v1.3.3"

    input:
    tuple val(meta), path(genbank)

    output:
    tuple val(meta), path("${prefix}.fna"), emit: fasta
    tuple val("${task.process}"), val("python"), eval("python3 --version 2>&1 | sed -n '1{s/Python //; p}'"), topic: versions, emit: versions_python


    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    extract_genome_fasta.py \\
        --genbank ${genbank} \\
        --output ${prefix}.fna
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    printf '>%s\\nACGTACGTACGTACGT\\n' "${prefix}" > ${prefix}.fna
    """
}
