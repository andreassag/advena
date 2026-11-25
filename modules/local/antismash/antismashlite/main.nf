/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ANTISMASH_ANTISMASHLITE process
    Detects secondary metabolite and bacteriocin BGCs on annotated GenBank genomes.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process ANTISMASH_ANTISMASHLITE {
    tag "${meta.id}"
    label "process_medium"

    conda "bioconda::antismash-lite=7.1.0"
    container "quay.io/biocontainers/antismash-lite:7.1.0--pyhdfd78af_0"

    input:
    tuple val(meta), path(genbank)
    path db

    output:
    tuple val(meta), path("${prefix}/${prefix}.json"), emit: json
    tuple val(meta), path("${prefix}/*.region*.gbk"), emit: region_gbks, optional: true
    tuple val("${task.process}"), val("antismash"), eval("antismash --version 2>&1 | sed -n '1{s/antiSMASH //; p}' || echo 'unknown'"), topic: versions, emit: versions_antismash


    script:
    def args = task.ext.args ?: "--minimal --genefinding-tool none"
    prefix = task.ext.prefix ?: "${meta.id}"
    def db_arg = db ? "--databases ${db}" : ""
    """
    gzip -cdf ${genbank} > ${prefix}.gbk
    mkdir -p ${prefix}

    antismash \\
        ${args} \\
        ${db_arg} \\
        --cpus ${task.cpus} \\
        --output-dir ${prefix} \\
        --output-basename ${prefix} \\
        ${prefix}.gbk
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}
    printf '{"records": []}\\n' > ${prefix}/${prefix}.json
    """
}
