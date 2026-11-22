/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    HMMER_HMMSCAN process
    Runs hmmscan against a set of pressed HMM profiles and formats 23-column TSV output.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process HMMER_HMMSCAN {
    tag "${meta.id}"
    label "process_medium"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/b1/b171806a7b780537073c1e9ffe50885ba7917d8072b5941a125698df1253c4f7/data'
        : 'community.wave.seqera.io/library/hmmer:3.4--ab3efa7c16f46920' }"

    input:
    tuple val(meta), path(faa)
    path hmm_files

    output:
    tuple val(meta), path("${prefix}_hmmscan_all.tsv"), emit: results
    tuple val(meta), path("per_profile/${prefix}_*.tsv"), emit: per_profile_tsv
    tuple val(meta), path("per_profile/${prefix}_*.txt"), emit: per_profile_txt
    tuple val("${task.process}"), val("hmmer"), eval("hmmscan -h 2>/dev/null | sed '2!d;s/^# HMMER *//;s/ .*//' || true"), topic: versions, emit: versions_hmmer

    script:
    def args = task.ext.args ?: ""
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    HEADER="target_name\\taccession\\ttlen\\tquery_name\\tseq_accession\\tqlen\\tseq_E-value\\tseq_score\\tseq_bias\\t#_domain\\tof_domain\\tdomain_c-Evalue\\tdomain_i-Evalue\\tdomain_score\\tdomain_bias\\thmm_coord_from\\thmm_coord_to\\tali_coord_from\\tali_coord_to\\tenv_coord_from\\tenv_coord_to\\tenv_coord_acc\\tdescription_of_target"

    mkdir -p per_profile

    first=true
    for hmmfile in \$(ls *.hmm 2>/dev/null | sort); do
        hmmname=\$(basename "\$hmmfile" .hmm)

        if [ ! -s "${faa}" ]; then
            touch "per_profile/${prefix}_\${hmmname}.txt"
            touch "per_profile/${prefix}_\${hmmname}.tsv"
        else
            hmmscan \\
                ${args} \\
                --domtblout "per_profile/${prefix}_\${hmmname}.txt" \\
                --cpu ${task.cpus} \\
                "\$hmmfile" ${faa} > /dev/null

            sed '/^#/d' "per_profile/${prefix}_\${hmmname}.txt" | \\
            awk -F '[[:space:]]+' 'BEGIN{ OFS="\\t" } { descr=\$23; for(i=24;i<=NF;i++) descr=descr" "\$i; print \$1,\$2,\$3,\$4,\$5,\$6,\$7,\$8,\$9,\$10,\$11,\$12,\$13,\$14,\$15,\$16,\$17,\$18,\$19,\$20,\$21,\$22,descr }' \\
            > "per_profile/${prefix}_\${hmmname}.tsv"
        fi

        if [ "\$first" = true ]; then
            echo -e "\$HEADER" > ${prefix}_hmmscan_all.tsv
            cat "per_profile/${prefix}_\${hmmname}.tsv" >> ${prefix}_hmmscan_all.tsv
            first=false
        else
            cat "per_profile/${prefix}_\${hmmname}.tsv" >> ${prefix}_hmmscan_all.tsv
        fi
    done

    if [ "\$first" = true ]; then
        echo -e "\$HEADER" > ${prefix}_hmmscan_all.tsv
    fi
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    HEADER="target_name\\taccession\\ttlen\\tquery_name\\tseq_accession\\tqlen\\tseq_E-value\\tseq_score\\tseq_bias\\t#_domain\\tof_domain\\tdomain_c-Evalue\\tdomain_i-Evalue\\tdomain_score\\tdomain_bias\\thmm_coord_from\\thmm_coord_to\\tali_coord_from\\tali_coord_to\\tenv_coord_from\\tenv_coord_to\\tenv_coord_acc\\tdescription_of_target"
    echo -e "\$HEADER" > ${prefix}_hmmscan_all.tsv

    mkdir -p per_profile
    found=false
    for hmmfile in \$(ls *.hmm 2>/dev/null | sort); do
        hmmname=\$(basename "\$hmmfile" .hmm)
        touch "per_profile/${prefix}_\${hmmname}.tsv"
        touch "per_profile/${prefix}_\${hmmname}.txt"
        found=true
    done

    if [ "\$found" = false ]; then
        touch "per_profile/${prefix}_stub.tsv"
        touch "per_profile/${prefix}_stub.txt"
    fi
    """
}
