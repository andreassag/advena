/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RGI_MAIN process
    Runs CARD Resistance Gene Identifier (RGI) in protein mode (Strict + Perfect by default).
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process RGI_MAIN {
    tag "${meta.id}"
    label "process_low"

    conda "bioconda::rgi=6.0.3"
    container "quay.io/biocontainers/rgi:6.0.3--pyha8f3691_1"

    input:
    tuple val(meta), path(faa)
    path card_db

    output:
    tuple val(meta), path("${prefix}_card_rgi.txt"), emit: tsv
    tuple val(meta), path("${prefix}_card_rgi.json"), emit: json
    tuple val("${task.process}"), val("rgi"), eval("rgi main --version 2>&1 | sed -n '1p' || echo 'unknown'"), topic: versions, emit: versions_rgi


    script:
    def args = task.ext.args ?: ""
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    export HOME="\$PWD"
    gzip -cdf ${faa} > query.faa

    if [ ! -s query.faa ]; then
        printf 'ORF_ID\\tContig\\tStart\\tStop\\tOrientation\\tCut_Off\\tBest_Hit_ARO\\tBest_Identities\\tARO\\tDrug Class\\tResistance Mechanism\\tPercentage Length of Reference Sequence\\n' > ${prefix}_card_rgi.txt
        printf '{}' > ${prefix}_card_rgi.json
    else
        CARD_JSON=\$(find -L ${card_db} -maxdepth 2 -name "card.json" | sed -n '1p')
        if [ -z "\${CARD_JSON}" ]; then
            echo "ERROR: card.json not found under ${card_db}" >&2
            exit 1
        fi
        rgi load --card_json "\${CARD_JSON}" --local
        rgi main \\
            ${args} \\
            --input_sequence query.faa \\
            --output_file ${prefix}_card_rgi \\
            --input_type protein \\
            --alignment_tool DIAMOND \\
            --num_threads ${task.cpus} \\
            --local \\
            --clean
    fi
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    printf 'ORF_ID\\tContig\\tStart\\tStop\\tOrientation\\tCut_Off\\tBest_Hit_ARO\\tBest_Identities\\tARO\\tDrug Class\\tResistance Mechanism\\tPercentage Length of Reference Sequence\\n' > ${prefix}_card_rgi.txt
    printf '{}' > ${prefix}_card_rgi.json
    """
}
