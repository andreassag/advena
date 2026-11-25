/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CARD_DOWNLOAD process
    Downloads and caches the latest Comprehensive Antibiotic Resistance Database (card.json).
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process CARD_DOWNLOAD {
    tag "card_db"
    label "process_single"

    conda "conda-forge::python=3.12"
    container "quay.io/biocontainers/rgi:6.0.3--pyha8f3691_1"

    output:
    path "card_db", emit: db
    tuple val("${task.process}"), val("card"), eval("date +%Y-%m-%d || echo 'unknown'"), topic: versions, emit: versions_card


    script:
    """
    mkdir -p card_db
    python3 - <<'PYEOF'
import os
import shutil
import tarfile
import urllib.request

archive = "card_data.tar.bz2"
urllib.request.urlretrieve("https://card.mcmaster.ca/latest/data", archive)
extracted = False
with tarfile.open(archive, "r:bz2") as tar:
    for member in tar.getmembers():
        if os.path.basename(member.name) == "card.json" and member.isfile():
            src = tar.extractfile(member)
            if src is not None:
                with open("card_db/card.json", "wb") as dst:
                    shutil.copyfileobj(src, dst)
                extracted = True
                break
if not extracted:
    raise RuntimeError("card.json not found in CARD archive")
PYEOF
    rm -f card_data.tar.bz2
    """

    stub:
    """
    mkdir -p card_db
    printf '{"_version": "stub"}\\n' > card_db/card.json
    """
}
