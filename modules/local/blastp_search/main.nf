/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    BLASTP_SEARCH process
    Runs BLASTP to search for signature proteins against one database.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process BLASTP_SEARCH {
    tag "${meta.id}_${db_name}"
    label "process_medium"
    conda "${moduleDir}/../environment.yml"
    container "ghcr.io/andreassag/icescreen:v1.3.3"

    input:
    tuple val(meta), path(faa), val(db_name)
    path db

    output:
    tuple val(meta), val(db_name), path("${meta.id}_${db_name}.tsv"), emit: results
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def db_path = "${db}/blastdb/${db_name}"
    """
    # Run BLASTP
    blastp \\
        -out ${meta.id}_${db_name}.tsv \\
        -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen' \\
        -db ${db_path} \\
        -evalue ${params.blastp_evalue} \\
        -query ${faa} \\
        -task blastp \\
        -max_target_seqs ${params.blastp_max_target_seqs} \\
        -num_threads ${task.cpus} || true

    # Prepend column header to the tabular BLAST output
    blastpfields="qseqid\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tqlen\tslen"

    if [ -s ${meta.id}_${db_name}.tsv ]; then
        { echo -e "\$blastpfields"; cat ${meta.id}_${db_name}.tsv; } > tmp_blast_header.tsv
        mv tmp_blast_header.tsv ${meta.id}_${db_name}.tsv
    else
        echo -e "\$blastpfields" > ${meta.id}_${db_name}.tsv
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        blast: \$(blastp -version 2>&1 | head -1 | sed 's/blastp: //')
    END_VERSIONS
    """

    stub:
    """
    echo -e "qseqid\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tqlen\tslen" > ${meta.id}_${db_name}.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        blast: 2.16.0
    END_VERSIONS
    """
}
