/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    INTEGRONFINDER process
    Detects complete integrons, CALINs, and attC cassettes in genome assemblies.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process INTEGRONFINDER {
    tag "${meta.id}"
    label "process_medium"

    conda "bioconda::integron_finder=2.0.5"
    container "quay.io/biocontainers/integron_finder:2.0.5--pyhdfd78af_0"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("${prefix}/${prefix}.integrons"), emit: integrons
    tuple val(meta), path("${prefix}/${prefix}.summary"), emit: summary
    tuple val("${task.process}"), val("integronfinder"), eval("integron_finder --version 2>&1 | sed -n '1{s/integron_finder version //; p}' || echo 'unknown'"), topic: versions, emit: versions_integronfinder


    script:
    def args = task.ext.args ?: "--local-max --mute"
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    gzip -cdf ${fasta} > ${prefix}.fna

    integron_finder \\
        ${args} \\
        --cpu ${task.cpus} \\
        --outdir ${prefix}_raw \\
        ${prefix}.fna

    mkdir -p ${prefix}
    INT_FILE=\$(find ${prefix}_raw -name "*.integrons" | sed -n '1p')
    SUM_FILE=\$(find ${prefix}_raw -name "*.summary" | sed -n '1p')
    if [ -n "\${INT_FILE}" ] && [ -f "\${INT_FILE}" ]; then
        cp "\${INT_FILE}" ${prefix}/${prefix}.integrons
    else
        printf 'ID_integron\\tID_replicon\\telement\\tpos_beg\\tpos_end\\tstrand\\tevalue\\ttype_elt\\tannotation\\tmodel\\ttype\\tdefault\\tdistance_2attC\\tconsidered_topology\\n' > ${prefix}/${prefix}.integrons
    fi
    if [ -n "\${SUM_FILE}" ] && [ -f "\${SUM_FILE}" ]; then
        cp "\${SUM_FILE}" ${prefix}/${prefix}.summary
    else
        printf 'ID_replicon\\tCALIN\\tcomplete\\tIn0\\ttopology\\tsize\\n' > ${prefix}/${prefix}.summary
    fi
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}
    printf 'ID_integron\\tID_replicon\\tpos_beg\\tpos_end\\n' > ${prefix}/${prefix}.integrons
    printf 'ID_replicon\\tCALIN\\tcomplete\\tIn0\\n' > ${prefix}/${prefix}.summary
    """
}
