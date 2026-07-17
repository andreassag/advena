

process DOWNLOAD_DB {
    // If db_dir is specified, use it. If false/null, default to the normal work/ directory.
    storeDir params.db_dir ? "${params.db_dir}" : null

    output:
    path "icescreen_database", emit: db

    script:
    """
    # Download and extract in a single pipe to save space, or download then extract
    wget -qO- "${params.db_url}" | tar -xzf -

    mv database icescreen_database

    find icescreen_database -name "*.hmm" -exec hmmpress -f {} \\;
    """

    stub:
    """
    mkdir -p icescreen_database
    """
}
