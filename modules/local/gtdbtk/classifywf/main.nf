/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    GTDBTK_CLASSIFYWF process
    Performs GTDB-Tk taxonomic classification to verify Bacillota phylum assignment.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process GTDBTK_CLASSIFYWF {
    tag "${meta.id}"
    label "process_high"

    conda "bioconda::gtdbtk=2.4.0"
    container "quay.io/biocontainers/gtdbtk:2.4.0--pyhdfd78af_1"

    input:
    tuple val(meta), path(fasta)
    path db

    output:
    tuple val(meta), path("${prefix}/classify/${prefix}.bac120.summary.tsv"), emit: summary
    tuple val("${task.process}"), val("gtdbtk"), eval("gtdbtk --version 2>&1 | sed -n '1{s/gtdbtk: version //; s/ .*//; p}' || echo 'unknown'"), topic: versions, emit: versions_gtdbtk


    script:
    def args = task.ext.args ?: "--skip_ani_screen"
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    export HOME="\$PWD"
    GTDB_ROOT=\$(find -L ${db} -maxdepth 2 -type d -name "taxonomy" -exec dirname {} \\; | sed -n '1p')
    export GTDBTK_DATA_PATH="\${GTDB_ROOT:-${db}}"
    mkdir -p genome_dir ${prefix}/classify
    gzip -cdf ${fasta} > genome_dir/${meta.id}.fna

    gtdbtk classify_wf \\
        ${args} \\
        --genome_dir genome_dir \\
        --extension fna \\
        --out_dir ${prefix} \\
        --prefix ${prefix} \\
        --cpus ${task.cpus}

    if [ ! -f ${prefix}/classify/${prefix}.bac120.summary.tsv ]; then
        printf 'user_genome\\tclassification\\n%s\\tUnclassified\\n' "${meta.id}" > ${prefix}/classify/${prefix}.bac120.summary.tsv
    fi
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}/classify
    printf 'user_genome\\tclassification\\n%s\\td__Bacteria;p__Bacillota;c__Bacilli\\n' "${meta.id}" > ${prefix}/classify/${prefix}.bac120.summary.tsv
    """
}
