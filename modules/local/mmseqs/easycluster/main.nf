/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    MMSEQS_EASYCLUSTER process
    Clusters all pan-ME protein sequences across the cohort using MMseqs2 easy-cluster.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process MMSEQS_EASYCLUSTER {
    tag "cohort_mmseqs"
    label "process_medium"

    conda "bioconda::mmseqs2=17.b804f"
    container "quay.io/biocontainers/mmseqs2:17.b804f--hd6d6fdc_1"

    input:
    path me_faas, stageAs: "me_faas/*"

    output:
    path "pan_me_cluster.tsv", emit: clusters
    path "pan_me_rep_seq.fasta", emit: rep_seqs
    tuple val("${task.process}"), val("mmseqs2"), eval("mmseqs version 2>&1 | sed -n '1p' || echo 'unknown'"), topic: versions, emit: versions_mmseqs2


    script:
    def args = task.ext.args ?: "--min-seq-id 0.8 -c 0.8"
    """
    : > combined_me_cds.faa
    for f in me_faas/*; do
        [ -e "\$f" ] || continue
        fname="\${f##*/}"
        sample="\${fname%_ME_cds.faa.gz}"
        sample="\${sample%_ME_cds.faa}"
        gzip -cdf "\$f" | awk -v s="\${sample}" '
            /^>/ {
                sub(/^>/, "")
                print ">" s "|" \$0
                next
            }
            { print }
        ' >> combined_me_cds.faa
    done

    if [ ! -s combined_me_cds.faa ]; then
        touch pan_me_cluster.tsv pan_me_rep_seq.fasta
    else
        mkdir -p mmseqs_tmp
        mmseqs easy-cluster \\
            combined_me_cds.faa \\
            pan_me \\
            mmseqs_tmp \\
            --threads ${task.cpus} \\
            ${args}
        rm -rf mmseqs_tmp
    fi
    """

    stub:
    """
    touch pan_me_cluster.tsv pan_me_rep_seq.fasta
    """
}
