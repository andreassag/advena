/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    GTDBTK_DOWNLOAD process
    Downloads and caches the GTDB-Tk reference database (>80 GB) when
    --auto_download_large_dbs true is explicitly enabled.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process GTDBTK_DOWNLOAD {
    tag "gtdbtk_db"
    label "process_low"

    conda "bioconda::gtdbtk=2.4.0"
    container "quay.io/biocontainers/gtdbtk:2.4.0--pyhdfd78af_1"

    output:
    path "gtdbtk_db", emit: db
    tuple val("${task.process}"), val("gtdbtk"), eval("gtdbtk --version 2>&1 | sed -n '1{s/gtdbtk: version //; s/ .*//; p}' || echo 'unknown'"), topic: versions, emit: versions_gtdbtk


    script:
    """
    mkdir -p gtdbtk_db
    python3 - <<'PYEOF'
import tarfile
import urllib.request

url = "https://data.gtdb.ecogenomic.org/releases/latest/auxillary_files/gtdbtk_package/full_package/gtdbtk_data.tar.gz"
archive = "gtdbtk_data.tar.gz"
urllib.request.urlretrieve(url, archive)
with tarfile.open(archive, "r:gz") as tar:
    tar.extractall("gtdbtk_db")
PYEOF
    rm -f gtdbtk_data.tar.gz
    """

    stub:
    """
    mkdir -p gtdbtk_db/release_stub
    touch gtdbtk_db/release_stub/gtdbtk_db.txt
    """
}
