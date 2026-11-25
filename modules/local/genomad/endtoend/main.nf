/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    GENOMAD_ENDTOEND process
    Identifies plasmid contigs and provirus/prophage regions using geNomad.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process GENOMAD_ENDTOEND {
    tag "${meta.id}"
    label "process_medium"

    conda "bioconda::genomad=1.8.1"
    container "quay.io/biocontainers/genomad:1.8.1--pyhdfd78af_0"

    input:
    tuple val(meta), path(fasta)
    path db

    output:
    tuple val(meta), path("${prefix}/${prefix}_plasmid_summary.tsv"), emit: plasmid_summary
    tuple val(meta), path("${prefix}/${prefix}_virus_summary.tsv"), emit: virus_summary
    tuple val("${task.process}"), val("genomad"), eval("genomad --version 2>&1 | sed -n '1{s/geNomad, version //; p}' || echo 'unknown'"), topic: versions, emit: versions_genomad


    script:
    def args = task.ext.args ?: "--cleanup"
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    gzip -cdf ${fasta} > ${prefix}.fna

    GENOMAD_DB=\$(find -L ${db} -maxdepth 2 -type f -name "genomad_db" -exec dirname {} \\; | sed -n '1p')
    if [ -z "\${GENOMAD_DB}" ]; then
        GENOMAD_DB=\$(find -L ${db} -maxdepth 2 -type f -name "version.txt" -exec dirname {} \\; | sed -n '1p')
    fi
    GENOMAD_DB=\${GENOMAD_DB:-${db}}

    genomad end-to-end \\
        ${args} \\
        --threads ${task.cpus} \\
        ${prefix}.fna \\
        ${prefix}_raw \\
        "\${GENOMAD_DB}"

    mkdir -p ${prefix}
    cp ${prefix}_raw/${prefix}_summary/${prefix}_plasmid_summary.tsv ${prefix}/${prefix}_plasmid_summary.tsv
    cp ${prefix}_raw/${prefix}_summary/${prefix}_virus_summary.tsv ${prefix}/${prefix}_virus_summary.tsv
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}
    printf 'seq_name\\tlength\\ttopology\\tn_genes\\tgenetic_code\\tplasmid_score\\tfdr\\tn_hallmarks\\tmarker_enrichment\\tconjugation_genes\\tamr_genes\\n' > ${prefix}/${prefix}_plasmid_summary.tsv
    printf 'seq_name\\tlength\\ttopology\\tcoordinates\\tn_genes\\tgenetic_code\\tvirus_score\\tfdr\\tn_hallmarks\\tmarker_enrichment\\ttaxonomy\\n' > ${prefix}/${prefix}_virus_summary.tsv
    """
}
