/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    DEFENSEFINDER_RUN process
    Detects anti-phage defense and toxin-antitoxin systems on ordered element proteomes.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process DEFENSEFINDER_RUN {
    tag "${meta.id}"
    label "process_low"

    conda "bioconda::defense-finder=1.3.0"
    container "quay.io/biocontainers/defense-finder:1.3.0--pyhdfd78af_0"

    input:
    tuple val(meta), path(faa)
    path models_dir

    output:
    tuple val(meta), path("${prefix}_defense_finder_systems.tsv"), emit: systems
    tuple val(meta), path("${prefix}_defense_finder_genes.tsv"), emit: genes
    tuple val("${task.process}"), val("defensefinder"), eval("defense-finder version 2>&1 | sed -n '1p' || echo 'unknown'"), topic: versions, emit: versions_defensefinder


    script:
    def args = task.ext.args ?: "--preserve-raw"
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    gzip -cdf ${faa} > query.faa

    if [ ! -s query.faa ]; then
        printf 'sys_id\\ttype\\tsubtype\\tactivity\\tsys_beg\\tsys_end\\tprotein_in_syst\\tgenes_count\\tname_of_profiles_in_sys\\n' > ${prefix}_defense_finder_systems.tsv
        printf 'replicon\\thit_id\\tgene_name\\thit_pos\\tmodel_fqn\\tsys_id\\tsys_loci\\tlocus_num\\tsys_wholeness\\tsys_score\\tsys_occ\\thit_gene_ref\\thit_status\\thit_seq_len\\thit_i_eval\\thit_score\\thit_profile_cov\\thit_seq_cov\\thit_begin_match\\thit_end_match\\tcounterpart\\tused_in\\ttype\\tsubtype\\tactivity\\n' > ${prefix}_defense_finder_genes.tsv
    else
        mkdir -p df_out
        defense-finder run \\
            ${args} \\
            --models-dir ${models_dir} \\
            --workers ${task.cpus} \\
            --out-dir df_out \\
            query.faa

        cp df_out/*defense_finder_systems.tsv ${prefix}_defense_finder_systems.tsv
        cp df_out/*defense_finder_genes.tsv ${prefix}_defense_finder_genes.tsv
    fi
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    printf 'sys_id\\ttype\\tsubtype\\n' > ${prefix}_defense_finder_systems.tsv
    printf 'hit_id\\tsys_id\\tgene_name\\ttype\\tsubtype\\n' > ${prefix}_defense_finder_genes.tsv
    """
}
