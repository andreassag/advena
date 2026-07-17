/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FP_SCREENING process
    Screens for false positive signature proteins using HMM FP profiles.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process FP_SCREENING {
    tag "${meta.id}"
    label "process_low"
    conda "${moduleDir}/../environment.yml"
    container "ghcr.io/exterex/icescreen:v1.3.3"

    input:
    tuple val(meta), path(merged_sp), path(faa)
    path db

    output:
    tuple val(meta), path("${meta.id}_detected_SP_hmm_cleaned.tsv"), emit: sp_clean
    tuple val(meta), path("${meta.id}_detected_SP_source.faa"), emit: sp_faa
    tuple val(meta), path("${meta.id}_fp_all.tsv"), emit: fp_all
    tuple val(meta), path("${meta.id}_hmm_fp_hits.tsv"), emit: fp_hits
    tuple val(meta), path("fp_per_profile/${meta.id}_*.tsv"), emit: fp_per_profile_tsv
    tuple val(meta), path("fp_per_profile/${meta.id}_*.txt"), emit: fp_per_profile_txt
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def mode_file = "${db}/mode/${params.phylum}.yml"
    def fp_profiles_dir = "${db}/hmmdb/FP_profiles"
    """
    HEADER="target_name\\taccession\\ttlen\\tquery_name\\tseq_accession\\tqlen\\tseq_E-value\\tseq_score\\tseq_bias\\t#_domain\\tof_domain\\tdomain_c-Evalue\\tdomain_i-Evalue\\tdomain_score\\tdomain_bias\\thmm_coord_from\\thmm_coord_to\\tali_coord_from\\tali_coord_to\\tenv_coord_from\\tenv_coord_to\\tenv_coord_acc\\tdescription_of_target"

    # Extract SP protein sequences from the FAA file using Python
    python3 - << 'PYEOF'
import sys, csv
from pathlib import Path

merged_tsv = "${merged_sp}"
faa_in    = "${faa}"
faa_out   = "${meta.id}_detected_SP_source.faa"

protein_ids = set()
with open(merged_tsv) as f:
    reader = csv.reader(f, delimiter="\\t")
    next(reader, None)  # skip header
    for row in reader:
        if row and len(row) > 2:
            protein_ids.add(row[2])

if not protein_ids:
    Path(faa_out).touch()
    sys.exit(0)

with open(faa_in) as f_in, open(faa_out, "w") as f_out:
    write = False
    for line in f_in:
        if line.startswith(">"):
            write = line[1:].strip() in protein_ids
        if write:
            f_out.write(line)
PYEOF

    mkdir -p fp_per_profile

    # Run hmmscan against each FP profile
    first=true
    for hmmfile in ${fp_profiles_dir}/*.hmm; do
        hmmname=\$(basename "\$hmmfile" .hmm)

        nrow_faa=\$(wc -l < ${meta.id}_detected_SP_source.faa)
        if [ "\$nrow_faa" -eq 0 ]; then
            touch "fp_per_profile/${meta.id}_\${hmmname}.txt"
            touch "fp_per_profile/${meta.id}_\${hmmname}.tsv"
        else
            hmmscan --domtblout "fp_per_profile/${meta.id}_\${hmmname}.txt" \\
                --cpu ${task.cpus} \\
                "\$hmmfile" ${meta.id}_detected_SP_source.faa > /dev/null

            sed '/^#/d' "fp_per_profile/${meta.id}_\${hmmname}.txt" | \\
            awk -F '[[:space:]]+' 'BEGIN{ OFS="\\t" } { descr=\$23; for(i=24;i<=NF;i++) descr=descr" "\$i; print \$1,\$2,\$3,\$4,\$5,\$6,\$7,\$8,\$9,\$10,\$11,\$12,\$13,\$14,\$15,\$16,\$17,\$18,\$19,\$20,\$21,\$22,descr }' \\
            > "fp_per_profile/${meta.id}_\${hmmname}.tsv"
        fi

        if [ "\$first" = true ]; then
            echo -e "\$HEADER" > ${meta.id}_fp_all.tsv
            cat "fp_per_profile/${meta.id}_\${hmmname}.tsv" >> ${meta.id}_fp_all.tsv
            first=false
        else
            cat "fp_per_profile/${meta.id}_\${hmmname}.tsv" >> ${meta.id}_fp_all.tsv
        fi
    done

    if [ "\$first" = true ]; then
        echo -e "\$HEADER" > ${meta.id}_fp_all.tsv
    fi

    # Filter false positives
    process_hmmscan_fp \\
        --insp ${merged_sp} \\
        --infp ${meta.id}_fp_all.tsv \\
        -c ${mode_file} \\
        --outfp ${meta.id}_hmm_fp_hits.tsv \\
        --outfiltered ${meta.id}_detected_SP_hmm_cleaned.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hmmer: \$(hmmscan -h | grep '^# HMMER' | sed 's/# HMMER //' | sed 's/ .*//')
        python: \$(python3 --version 2>&1 | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}_detected_SP_hmm_cleaned.tsv
    touch ${meta.id}_detected_SP_source.faa
    touch ${meta.id}_fp_all.tsv
    touch ${meta.id}_hmm_fp_hits.tsv

    mkdir -p fp_per_profile
    for profile in FP_Pfam_a FP_Pfam_b FP_Pfam_c; do
        touch "fp_per_profile/${meta.id}_\${profile}.tsv"
        touch "fp_per_profile/${meta.id}_\${profile}.txt"
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hmmer: 3.4
        python: 3.11.0
    END_VERSIONS
    """
}
