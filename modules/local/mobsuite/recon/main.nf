/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    MOBSUITE_RECON process
    Classifies contigs into chromosome vs. plasmid and reconstructs plasmid clusters.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process MOBSUITE_RECON {
    tag "${meta.id}"
    label "process_medium"

    conda "bioconda::mob_suite=3.1.9"
    container "quay.io/biocontainers/mob_suite:3.1.9--pyhdfd78af_1"

    input:
    tuple val(meta), path(fasta)
    path db

    output:
    tuple val(meta), path("${prefix}/contig_report.txt"), emit: contig_report
    tuple val(meta), path("${prefix}/mobtyper_results.txt"), emit: mobtyper_results, optional: true
    tuple val(meta), path("${prefix}/plasmid_*.fasta"), emit: plasmids, optional: true
    tuple val("${task.process}"), val("mobsuite"), eval("mob_recon --version 2>&1 | sed -n '1{s/mob_recon //; p}' || echo 'unknown'"), topic: versions, emit: versions_mobsuite


    script:
    def args = task.ext.args ?: "--force"
    prefix = task.ext.prefix ?: "${meta.id}"
    def db_arg = db ? "--database_directory ${db}" : ""
    """
    gzip -cdf ${fasta} > ${prefix}.fna

    mob_recon \\
        ${args} \\
        ${db_arg} \\
        --num_threads ${task.cpus} \\
        --infile ${prefix}.fna \\
        --outdir ${prefix}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}
    printf 'sample_id\\tmolecule_type\\tprimary_cluster_id\\tcontig_id\\n%s\\tchromosome\\t-\\tNC_000964.3\\n' "${meta.id}" > ${prefix}/contig_report.txt
    """
}
