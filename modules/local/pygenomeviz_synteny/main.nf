/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PYGENOMEVIZ_SYNTENY process
    Renders SVG and interactive HTML synteny alignments per multi-member ICE/IME cluster.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process PYGENOMEVIZ_SYNTENY {
    tag "pygenomeviz_synteny"
    label "process_low"

    conda "bioconda::pygenomeviz=1.5.0"
    container "quay.io/biocontainers/pygenomeviz:1.5.0--pyhdfd78af_0"

    input:
    path clusters_tsv
    path gbk_slices, stageAs: "gbks/*"
    path cargo_profiles, stageAs: "cargo/*"

    output:
    path "synteny_plots/*_synteny.svg", emit: svgs, optional: true
    path "synteny_plots/*_synteny.html", emit: htmls, optional: true
    path "synteny_plots/synteny_manifest.txt", emit: manifest
    tuple val("${task.process}"), val("pygenomeviz"), eval("python3 -c 'import pygenomeviz; print(pygenomeviz.__version__)' 2>/dev/null || echo 'standalone-svg'"), topic: versions, emit: versions_pygenomeviz


    script:
    def args = task.ext.args ?: ""
    def gbk_arg = gbk_slices ? "--gbks ${gbk_slices}" : ""
    def cargo_arg = cargo_profiles ? "--cargo-profiles ${cargo_profiles}" : ""
    """
    export HOME="\$PWD"
    export MPLBACKEND=Agg
    export MPLCONFIGDIR="\$PWD/.mplconfig"
    export XDG_CACHE_HOME="\$PWD/.cache"
    mkdir -p synteny_plots "\${MPLCONFIGDIR}" "\${XDG_CACHE_HOME}"
    plot_me_synteny.py \\
        --clusters ${clusters_tsv} \\
        ${gbk_arg} \\
        ${cargo_arg} \\
        ${args} \\
        --outdir synteny_plots
    ls -1 synteny_plots/*_synteny.svg 2>/dev/null > synteny_plots/synteny_manifest.txt || touch synteny_plots/synteny_manifest.txt
    """

    stub:
    """
    mkdir -p synteny_plots
    printf '<svg xmlns="http://www.w3.org/2000/svg"></svg>\\n' > synteny_plots/CLUSTER_001_synteny.svg
    printf '<!DOCTYPE html><html><body>stub</body></html>\\n' > synteny_plots/CLUSTER_001_synteny.html
    printf 'synteny_plots/CLUSTER_001_synteny.svg\\n' > synteny_plots/synteny_manifest.txt
    """
}
