/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    EXTRACT_ME_SEQUENCES process
    Slices detected ICE/IME intervals into BED, FNA, CDS FAA/GFF3, manifest, and per-ME GBK.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process EXTRACT_ME_SEQUENCES {
    tag "${meta.id}"
    label "process_single"

    conda "${moduleDir}/../environment.yml"
    container "ghcr.io/andreassag/icescreen:v1.3.3"

    input:
    tuple val(meta), path(detected_me), path(sp_with_me), path(annotated_gb)

    output:
    tuple val(meta), path("${prefix}_ME_regions.bed"), emit: bed
    tuple val(meta), path("${prefix}_ME_regions.fna.gz"), emit: fna
    tuple val(meta), path("${prefix}_cargo_query.fna.gz"), emit: cargo_fna
    tuple val(meta), path("${prefix}_ME_cds.faa.gz"), emit: faa
    tuple val(meta), path("${prefix}_ME_cds.gff3.gz"), emit: gff3
    tuple val(meta), path("${prefix}_ME_cds_manifest.tsv"), emit: manifest
    tuple val(meta), path("${prefix}_me_count.txt"), emit: me_count
    tuple val(meta), path("slices/*.gbk.gz"), emit: gbk_slices, optional: true
    tuple val("${task.process}"), val("python"), eval("python3 --version 2>&1 | sed -n '1{s/Python //; p}'"), topic: versions, emit: versions_python


    script:
    def args = task.ext.args ?: ""
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    extract_me_sequences.py \\
        --sample ${meta.id} \\
        --prefix ${prefix} \\
        --detected-me ${detected_me} \\
        --sp-with-me ${sp_with_me} \\
        --genbank ${annotated_gb} \\
        --outdir . \\
        ${args}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p slices
    printf 'NC_000964.3\\t504\\t12229\\tID_1\\t0\\t+\\tComplete ICE (I+C+R+V)\\n' > ${prefix}_ME_regions.bed
    printf '>%s__ID_1\\nACGTACGTACGTACGT\\n' "${meta.id}" | gzip > ${prefix}_ME_regions.fna.gz
    cp ${prefix}_ME_regions.fna.gz ${prefix}_cargo_query.fna.gz
    printf '>BSU_04800\\nMKKLLVVAA\\n' | gzip > ${prefix}_ME_cds.faa.gz
    printf '##gff-version 3\\n%s__ID_1\\tadvena\\tCDS\\t1\\t100\\t.\\t+\\t0\\tID=BSU_04800;locus_tag=BSU_04800\\n' "${meta.id}" | gzip > ${prefix}_ME_cds.gff3.gz
    printf 'sample\\tme_id\\tme_category\\tcontig\\tlocus_tag\\tprotein_id\\torig_start\\torig_end\\tstrand\\tis_signature_protein\\tsp_type\\tproduct\\n%s\\tID_1\\tComplete ICE\\tNC_000964.3\\tBSU_04800\\tNP_1\\t505\\t1611\\t-\\tTRUE\\tTyrosine integrase\\tintegrase\\n' "${meta.id}" > ${prefix}_ME_cds_manifest.tsv
    printf '1\\n' > ${prefix}_me_count.txt
    printf 'LOCUS       %s__ID_1              16 bp    DNA     linear   BCT 01-JAN-2026\\nFEATURES             Location/Qualifiers\\nORIGIN\\n//\\n' "${meta.id}" | gzip > slices/${meta.id}__ID_1.gbk.gz
    """
}
