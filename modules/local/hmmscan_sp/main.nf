/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    HMMSCAN_SP process
    Runs hmmscan against all signature protein HMM profiles.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process HMMSCAN_SP {
    tag "${meta.id}"
    label "process_medium"
    conda "${moduleDir}/../environment.yml" // TODO: Might not need the conda and container environment anymore for this process
    container "ghcr.io/andreassag/icescreen:v1.3.3"

    input:
    tuple val(meta), path(faa)
    path db

    output:
    tuple val(meta), path("${meta.id}_hmmscan_all.tsv"), emit: results
    tuple val(meta), path("per_profile/${meta.id}_*.tsv"), emit: per_profile_tsv
    tuple val(meta), path("per_profile/${meta.id}_*.txt"), emit: per_profile_txt
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def sp_profiles_dir = "${db}/hmmdb/SP_profiles"
    """
    HEADER="target_name\\taccession\\ttlen\\tquery_name\\tseq_accession\\tqlen\\tseq_E-value\\tseq_score\\tseq_bias\\t#_domain\\tof_domain\\tdomain_c-Evalue\\tdomain_i-Evalue\\tdomain_score\\tdomain_bias\\thmm_coord_from\\thmm_coord_to\\tali_coord_from\\tali_coord_to\\tenv_coord_from\\tenv_coord_to\\tenv_coord_acc\\tdescription_of_target"

    mkdir -p per_profile

    first=true
    for hmmfile in ${sp_profiles_dir}/*.hmm; do
        hmmname=\$(basename "\$hmmfile" .hmm)

        # Run hmmscan, store domtblout as per-profile .txt
        hmmscan --domtblout "per_profile/${meta.id}_\${hmmname}.txt" \\
            --cpu ${task.cpus} \\
            "\$hmmfile" ${faa} > /dev/null

        # Format to per-profile TSV
        sed '/^#/d' "per_profile/${meta.id}_\${hmmname}.txt" | \\
        awk -F '[[:space:]]+' 'BEGIN{ OFS="\\t"; } { descr=\$23; for(i=24; i<=NF; i++) descr=descr" "\$i; print \$1,\$2,\$3,\$4,\$5,\$6,\$7,\$8,\$9,\$10,\$11,\$12,\$13,\$14,\$15,\$16,\$17,\$18,\$19,\$20,\$21,\$22,descr; }' > "per_profile/${meta.id}_\${hmmname}.tsv"

        if [ "\$first" = true ]; then
            echo -e "\$HEADER" > ${meta.id}_hmmscan_all.tsv
            cat "per_profile/${meta.id}_\${hmmname}.tsv" >> ${meta.id}_hmmscan_all.tsv
            first=false
        else
            cat "per_profile/${meta.id}_\${hmmname}.tsv" >> ${meta.id}_hmmscan_all.tsv
        fi
    done

    # If no profiles found, create empty file with header
    if [ "\$first" = true ]; then
        echo -e "\$HEADER" > ${meta.id}_hmmscan_all.tsv
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hmmer: \$(hmmscan -h | grep '^# HMMER' | sed 's/# HMMER //' | sed 's/ .*//')
    END_VERSIONS
    """

    stub:
    """
    echo -e "target_name\\taccession\\ttlen\\tquery_name\\tseq_accession\\tqlen\\tseq_E-value\\tseq_score\\tseq_bias\\t#_domain\\tof_domain\\tdomain_c-Evalue\\tdomain_i-Evalue\\tdomain_score\\tdomain_bias\\thmm_coord_from\\thmm_coord_to\\tali_coord_from\\tali_coord_to\\tenv_coord_from\\tenv_coord_to\\tenv_coord_acc\\tdescription_of_target" > ${meta.id}_hmmscan_all.tsv

    mkdir -p per_profile
    for profile in Relaxase_PHA_IME_A1 Couplingprot_T4SS_t4cp1 Couplingprot2_T4SS_t4cp2 VirB4_T4SS_virb4 IntTyr_Phage_integrase IntSer_Recombinase DDE_UPF0236; do
        touch "per_profile/${meta.id}_\${profile}.tsv"
        touch "per_profile/${meta.id}_\${profile}.txt"
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hmmer: 3.4
    END_VERSIONS
    """
}
