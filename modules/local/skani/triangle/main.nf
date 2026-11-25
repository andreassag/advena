/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SKANI_TRIANGLE process
    Computes all-vs-all pairwise nucleotide ANI and aligned fraction across ICE/IME slices.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process SKANI_TRIANGLE {
    tag "cohort_skani"
    label "process_medium"

    conda "bioconda::skani=0.2.2"
    container "quay.io/biocontainers/skani:0.2.2--ha6fb395_2"

    input:
    path me_fnas, stageAs: "me_fnas/*"

    output:
    path "skani_triangle.tsv", emit: matrix
    tuple val("${task.process}"), val("skani"), eval("skani --version 2>&1 | sed -n '1{s/skani //; p}' || echo 'unknown'"), topic: versions, emit: versions_skani


    script:
    def args = task.ext.args ?: "--small-genomes -m 200 --min-af 50 --sparse"
    """
    mkdir -p split_elements
    for f in me_fnas/*; do
        [ -e "\$f" ] || continue
        gzip -cdf "\$f" | awk '
            /^>/ {
                if (out != "") close(out)
                hdr = substr(\$1, 2)
                gsub(/\\//, "_", hdr)
                out = "split_elements/" hdr ".fna"
                print ">" hdr > out
                next
            }
            out != "" && NF > 0 { print >> out }
            END { if (out != "") close(out) }
        '
    done

    N_ELEMS=\$(find split_elements -name "*.fna" -size +0c | wc -l)
    if [ "\${N_ELEMS}" -lt 2 ]; then
        printf 'Ref_file\\tQuery_file\\tANI\\tAlign_fraction_ref\\tAlign_fraction_query\\tRef_name\\tQuery_name\\n' > skani_triangle.tsv
    else
        skani triangle \\
            ${args} \\
            -t ${task.cpus} \\
            split_elements/*.fna \\
            -o skani_triangle.tsv
    fi
    """

    stub:
    """
    printf 'Ref_file\\tQuery_file\\tANI\\tAlign_fraction_ref\\tAlign_fraction_query\\tRef_name\\tQuery_name\\n' > skani_triangle.tsv
    """
}
