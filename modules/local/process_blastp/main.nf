/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PROCESS_BLASTP process
    Prepends ICEscreen tabular headers, filters BLASTP results, and selects best hits.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process PROCESS_BLASTP {
    tag "${meta.id}"
    label "process_low"

    conda "${moduleDir}/../environment.yml"
    container "ghcr.io/andreassag/icescreen:v1.3.3"

    input:
    tuple val(meta), path(blast_results)
    path mode_file
    path ice_finder_db

    output:
    tuple val(meta), path("${prefix}_blast_SP.tsv"), emit: best_hits
    tuple val(meta), path("${prefix}_IntTyr_filtered.tsv"), emit: intyr_hits
    tuple val(meta), path("blastp_output/${prefix}_*.tsv"), emit: blast_raw
    tuple val(meta), path("${prefix}_*_unfiltered.tsv"), emit: blast_unfiltered
    tuple val(meta), path("${prefix}_*_filtered.tsv"), emit: blast_filtered
    tuple val(meta), path("${prefix}_*_filtered_best.tsv"), emit: blast_filtered_best
    tuple val("${task.process}"), val("icescreen"), eval("python3 -c 'import importlib.metadata; print(importlib.metadata.version(\"icescreen\"))' 2>/dev/null || echo 'unknown'"), topic: versions, emit: versions_icescreen
    tuple val("${task.process}"), val("pandas"), eval("python3 -c 'import pandas; print(pandas.__version__)' 2>/dev/null || echo 'unknown'"), topic: versions, emit: versions_pandas

    script:
    def args = task.ext.args ?: ""
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    HEADER="qseqid\\tsseqid\\tpident\\tlength\\tmismatch\\tgapopen\\tqstart\\tqend\\tsstart\\tsend\\tevalue\\tbitscore\\tqlen\\tslen"
    mkdir -p blastp_output

    for blast_file in \$(ls *.tsv | sort); do
        db_name=\$(basename "\$blast_file" .tsv | sed "s/^${meta.id}_//")
        raw_out="blastp_output/${prefix}_\${db_name}.tsv"

        # Prepend 14-column header expected by ICEscreen if raw BLAST output is headerless
        if [ -s "\$blast_file" ] && head -n 1 "\$blast_file" | grep -q "^qseqid\t"; then
            cp "\$blast_file" "\$raw_out"
        elif [ -s "\$blast_file" ]; then
            { printf "%b\\n" "\$HEADER"; cat "\$blast_file"; } > "\$raw_out"
        else
            printf "%b\\n" "\$HEADER" > "\$raw_out"
        fi

        process_blastp_results \\
            ${args} \\
            -i "\$raw_out" \\
            -c ${mode_file} \\
            -d ${ice_finder_db} \\
            --outall "${prefix}_\${db_name}_unfiltered.tsv" \\
            --outfiltered "${prefix}_\${db_name}_filtered.tsv" \\
            --outbest "${prefix}_\${db_name}_filtered_best.tsv"
    done

    # Ensure IntTyr filtered file exists (may be absent if no IntTyr hits)
    if [ ! -f ${prefix}_IntTyr_filtered.tsv ]; then
        printf "%b\\n" "\$HEADER" > ${prefix}_IntTyr_filtered.tsv
    fi

    # Gather all best hits into one file
    first=true
    for best_file in \$(ls ${prefix}_*_filtered_best.tsv | sort); do
        if [ "\$first" = true ]; then
            cat "\$best_file" > all_best_combined.tsv
            first=false
        else
            tail -n+2 "\$best_file" >> all_best_combined.tsv
        fi
    done

    # Sort by CDS number and retain only the best hit per locus
    if [ -f all_best_combined.tsv ]; then
        tail -n+2 all_best_combined.tsv | sort -t\$'\\t' -k2n | cat <(head -1 all_best_combined.tsv) - > all_best_sorted.tsv
        retain_only_best_blast_hit_for_each_locus_tag \\
            -i all_best_sorted.tsv \\
            -o ${prefix}_blast_SP.tsv
    else
        touch ${prefix}_blast_SP.tsv
    fi
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    HEADER="qseqid\\tsseqid\\tpident\\tlength\\tmismatch\\tgapopen\\tqstart\\tqend\\tsstart\\tsend\\tevalue\\tbitscore\\tqlen\\tslen"
    mkdir -p blastp_output
    touch ${prefix}_blast_SP.tsv
    printf "%b\\n" "\$HEADER" > ${prefix}_IntTyr_filtered.tsv
    for db in Relaxase Couplingprot Couplingprot2 VirB4 IntTyr IntSer DDE; do
        printf "%b\\n" "\$HEADER" > "blastp_output/${prefix}_\${db}.tsv"
        touch "${prefix}_\${db}_unfiltered.tsv"
        touch "${prefix}_\${db}_filtered.tsv"
        touch "${prefix}_\${db}_filtered_best.tsv"
    done
    """
}
