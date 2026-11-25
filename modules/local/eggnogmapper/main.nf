/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    EGGNOGMAPPER process
    Runs eggNOG-mapper in DIAMOND mode to assign hierarchical OGs, COGs, KEGG, and PFAMs.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process EGGNOGMAPPER {
    tag "${meta.id}"
    label "process_medium"

    conda "bioconda::eggnog-mapper=2.1.12"
    container "quay.io/biocontainers/eggnog-mapper:2.1.12--pyhdfd78af_0"

    input:
    tuple val(meta), path(faa)
    path db

    output:
    tuple val(meta), path("${prefix}.emapper.annotations"), emit: annotations
    tuple val(meta), path("${prefix}.emapper.seed_orthologs"), emit: orthologs, optional: true
    tuple val("${task.process}"), val("eggnogmapper"), eval("emapper.py --version 2>&1 | sed -n '1{s/emapper-//;p}' || echo 'unknown'"), topic: versions, emit: versions_eggnogmapper


    script:
    def args = task.ext.args ?: ""
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    gzip -cdf ${faa} > query.faa

    if [ ! -s query.faa ]; then
        printf '#query\\tseed_ortholog\\tevalue\\tscore\\teggNOG_OGs\\tmax_annot_lvl\\tCOG_category\\tDescription\\tPreferred_name\\tGOs\\tEC\\tKEGG_ko\\tKEGG_Pathway\\tKEGG_Module\\tKEGG_Reaction\\tKEGG_rclass\\tBRITE\\tKEGG_TC\\tCAZy\\tBiGG_Reaction\\tPFAMs\\n' > ${prefix}.emapper.annotations
    else
        DB_DIR="${db}"
        if [ ! -f "\${DB_DIR}/eggnog_proteins.dmnd" ] && [ -f "\${DB_DIR}/eggnog_proteins.faa" ]; then
            cp -rL "\${DB_DIR}" eggnog_db_work
            DB_DIR="eggnog_db_work"
            diamond makedb --in "\${DB_DIR}/eggnog_proteins.faa" -d "\${DB_DIR}/eggnog_proteins" >/dev/null
        fi

        emapper.py \\
            ${args} \\
            -m diamond \\
            -i query.faa \\
            --data_dir "\${DB_DIR}" \\
            --cpu ${task.cpus} \\
            -o ${prefix}
    fi
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    printf '#query\\tseed_ortholog\\tevalue\\tscore\\teggNOG_OGs\\tmax_annot_lvl\\tCOG_category\\tDescription\\tPreferred_name\\tGOs\\tEC\\tKEGG_ko\\tKEGG_Pathway\\tKEGG_Module\\tKEGG_Reaction\\tKEGG_rclass\\tBRITE\\tKEGG_TC\\tCAZy\\tBiGG_Reaction\\tPFAMs\\n' > ${prefix}.emapper.annotations
    touch ${prefix}.emapper.seed_orthologs
    """
}
