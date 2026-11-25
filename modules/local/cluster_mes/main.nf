/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CLUSTER_MES process
    Assigns cohort ICE/IME cluster IDs using skani ANI/AF and MMseqs2 protein Jaccard.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process CLUSTER_MES {
    tag "cluster_mes"
    label "process_single"

    conda "${moduleDir}/../environment.yml"
    container "ghcr.io/andreassag/icescreen:v1.3.3"

    input:
    path beds, stageAs: "beds/*"
    path manifests, stageAs: "manifests/*"
    path skani_tsv
    path mmseqs_tsv

    output:
    path "me_clusters.tsv", emit: clusters_tsv
    path "cluster_members.json", emit: clusters_json
    tuple val("${task.process}"), val("python"), eval("python3 --version 2>&1 | sed -n '1{s/Python //; p}'"), topic: versions, emit: versions_python


    script:
    def args = task.ext.args ?: ""
    def beds_arg = beds ? "--beds ${beds}" : ""
    def man_arg = manifests ? "--manifests ${manifests}" : ""
    def skani_arg = skani_tsv ? "--skani ${skani_tsv}" : ""
    def mmseqs_arg = mmseqs_tsv ? "--mmseqs ${mmseqs_tsv}" : ""
    """
    cluster_mes.py \\
        ${beds_arg} \\
        ${man_arg} \\
        ${skani_arg} \\
        ${mmseqs_arg} \\
        ${args} \\
        --out-tsv me_clusters.tsv \\
        --out-json cluster_members.json
    """

    stub:
    """
    printf 'me_cluster_id\\tcluster_size\\telement_key\\tsample\\tme_id\\tme_category\\tlength_bp\\tn_cds\\tclustering_evidence\\n' > me_clusters.tsv
    printf '{}\\n' > cluster_members.json
    """
}
