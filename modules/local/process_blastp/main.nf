/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PROCESS_BLASTP process
    Filters BLASTP results and selects best hits per locus for all databases.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process PROCESS_BLASTP {
    tag "${meta.id}"
    label "process_low"
    conda "${moduleDir}/../environment.yml"
    container "ghcr.io/andreassag/icescreen:v1.3.3"

    input:
    tuple val(meta), val(db_names), path(blast_results)
    path db

    output:
    tuple val(meta), path("${meta.id}_blast_SP.tsv"), emit: best_hits
    tuple val(meta), path("${meta.id}_IntTyr_filtered.tsv"), emit: intyr_hits
    tuple val(meta), path("${meta.id}_*_unfiltered.tsv"), emit: blast_unfiltered
    tuple val(meta), path("${meta.id}_*_filtered.tsv"), emit: blast_filtered
    tuple val(meta), path("${meta.id}_*_filtered_best.tsv"), emit: blast_filtered_best
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def mode_file = "${db}/mode/${params.phylum}.yml"
    def ice_finder_db = "${db}/blastdb/ICE_Finder.db"
    """
    # Process each BLAST result file through filtering
    for blast_file in ${blast_results}; do
        db_name=\$(basename "\$blast_file" .tsv | sed "s/^${meta.id}_//")

        process_blastp_results \\
            -i "\$blast_file" \\
            -c ${mode_file} \\
            -d ${ice_finder_db} \\
            --outall "${meta.id}_\${db_name}_unfiltered.tsv" \\
            --outfiltered "${meta.id}_\${db_name}_filtered.tsv" \\
            --outbest "${meta.id}_\${db_name}_filtered_best.tsv"
    done

    # Ensure IntTyr filtered file exists (may be absent if no IntTyr hits)
    if [ ! -f ${meta.id}_IntTyr_filtered.tsv ]; then
        echo -e "qseqid\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tqlen\tslen" > ${meta.id}_IntTyr_filtered.tsv
    fi

    # Gather all best hits into one file
    first=true
    for best_file in ${meta.id}_*_filtered_best.tsv; do
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
            -o ${meta.id}_blast_SP.tsv
    else
        touch ${meta.id}_blast_SP.tsv
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version 2>&1 | sed 's/Python //')
        pandas: \$(python3 -c "import pandas; print(pandas.__version__)")
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}_blast_SP.tsv
    touch ${meta.id}_IntTyr_filtered.tsv
    for db in Relaxase Couplingprot Couplingprot2 VirB4 IntTyr IntSer DDE; do
        touch "${meta.id}_\${db}_unfiltered.tsv"
        touch "${meta.id}_\${db}_filtered.tsv"
        touch "${meta.id}_\${db}_filtered_best.tsv"
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: 3.11.0
        pandas: 2.2.0
    END_VERSIONS
    """
}
