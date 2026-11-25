/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    EXTRACT_SP_FAA process
    Extracts candidate signature protein sequences from a proteome FASTA using merged SP TSV.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process EXTRACT_SP_FAA {
    tag "${meta.id}"
    label "process_single"

    conda "${moduleDir}/../environment.yml"
    container "ghcr.io/andreassag/icescreen:v1.3.3"

    input:
    tuple val(meta), path(merged_sp), path(faa)

    output:
    tuple val(meta), path("${prefix}_detected_SP_source.faa"), emit: sp_faa
    tuple val("${task.process}"), val("python"), eval("python3 --version 2>&1 | sed -n '1{s/Python //; p}'"), topic: versions, emit: versions_python

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    python3 - << 'PYEOF'
import csv
import sys
from pathlib import Path

merged_tsv = "${merged_sp}"
faa_in = "${faa}"
faa_out = "${prefix}_detected_SP_source.faa"

protein_ids = set()
with open(merged_tsv, encoding="utf-8") as f:
    reader = csv.reader(f, delimiter="\\t")
    next(reader, None)  # skip header
    for row in reader:
        if row and len(row) > 2:
            protein_ids.add(row[2].strip())

if not protein_ids:
    Path(faa_out).touch()
    sys.exit(0)

with open(faa_in, encoding="utf-8") as f_in, open(faa_out, "w", encoding="utf-8") as f_out:
    write = False
    for line in f_in:
        if line.startswith(">"):
            header = line[1:].strip()
            cds_num = header.split("&", 1)[0]
            write = cds_num in protein_ids or header in protein_ids
        if write:
            f_out.write(line)
PYEOF
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_detected_SP_source.faa
    """
}
