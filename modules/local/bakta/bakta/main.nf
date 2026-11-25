/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    BAKTA_BAKTA process
    Annotates raw bacterial FASTA genome assemblies into GenBank, GFF3, and FAA outputs.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process BAKTA_BAKTA {
    tag "${meta.id}"
    label "process_medium"

    conda "bioconda::bakta=1.10.4"
    container "quay.io/biocontainers/bakta:1.10.4--pyhdfd78af_0"

    input:
    tuple val(meta), path(fasta)
    path db

    output:
    tuple val(meta), path("${prefix}/${prefix}.gbff"), emit: gbff
    tuple val(meta), path("${prefix}/${prefix}.gff3"), emit: gff3
    tuple val(meta), path("${prefix}/${prefix}.faa"), emit: faa
    tuple val(meta), path("${prefix}/${prefix}.txt"), emit: txt
    tuple val(meta), path("${prefix}/${prefix}.tsv"), emit: tsv
    tuple val("${task.process}"), val("bakta"), eval("bakta --version 2>&1 | sed -n '1{s/bakta //; p}' || echo 'unknown'"), topic: versions, emit: versions_bakta


    script:
    def args = task.ext.args ?: ""
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    DB_PATH=\$(find -L ${db} -maxdepth 2 -name "bakta.db" -exec dirname {} \\; | sed -n '1p')
    DB_PATH=\${DB_PATH:-${db}}

    bakta \\
        ${args} \\
        --db "\${DB_PATH}" \\
        --threads ${task.cpus} \\
        --prefix ${prefix} \\
        --output ${prefix} \\
        ${fasta}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}
    printf 'LOCUS       %s                   100 bp    DNA     linear   BCT 01-JAN-2026\\nFEATURES             Location/Qualifiers\\nORIGIN\\n//\\n' "${prefix}" > ${prefix}/${prefix}.gbff
    printf '##gff-version 3\\n' > ${prefix}/${prefix}.gff3
    touch ${prefix}/${prefix}.faa
    touch ${prefix}/${prefix}.txt
    touch ${prefix}/${prefix}.tsv
    """
}
