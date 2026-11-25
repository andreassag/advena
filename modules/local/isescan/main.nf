/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ISESCAN process
    Detects insertion sequence (IS) elements and families across genome assemblies.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process ISESCAN {
    tag "${meta.id}"
    label "process_medium"

    conda "bioconda::isescan=1.7.2.3"
    container "quay.io/biocontainers/isescan:1.7.2.3--h7b50bb2_3"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("${prefix}/${prefix}_isescan.tsv"), emit: results
    tuple val("${task.process}"), val("isescan"), eval("isescan.py --version 2>&1 | sed -n '1{s/ISEScan //; p}' || echo 'unknown'"), topic: versions, emit: versions_isescan


    script:
    def args = task.ext.args ?: ""
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    gzip -cdf ${fasta} > ${prefix}.fna

    isescan.py \\
        ${args} \\
        --seqfile ${prefix}.fna \\
        --output ${prefix}_raw \\
        --nthread ${task.cpus}

    mkdir -p ${prefix}
    TSV=\$(find ${prefix}_raw -name "*.tsv" | sed -n '1p')
    if [ -n "\${TSV}" ] && [ -s "\${TSV}" ]; then
        cp "\${TSV}" ${prefix}/${prefix}_isescan.tsv
    else
        printf 'seqID\\tfamily\\tcluster\\tisBegin\\tisEnd\\tisLen\\tncopy4is\\tstart1\\tend1\\tstart2\\tend2\\tscore\\tirId\\tirLen\\tnGaps\\torfBegin\\torfEnd\\tstrand\\torfLen\\tE-value\\tE-value4copy\\ttype\\tov\\ttir\\n' > ${prefix}/${prefix}_isescan.tsv
    fi
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}
    printf 'seqID\\tfamily\\tisBegin\\tisEnd\\n' > ${prefix}/${prefix}_isescan.tsv
    """
}
