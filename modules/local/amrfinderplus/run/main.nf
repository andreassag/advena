/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    AMRFINDERPLUS_RUN process
    Runs NCBI AMRFinderPlus (--plus) on extracted ME proteins, nucleotides, and GFF3.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process AMRFINDERPLUS_RUN {
    tag "${meta.id}"
    label "process_low"

    conda "bioconda::ncbi-amrfinderplus=4.0.19"
    container "quay.io/biocontainers/ncbi-amrfinderplus:4.0.19--h6e70893_0"

    input:
    tuple val(meta), path(faa), path(fna), path(gff)
    path db

    output:
    tuple val(meta), path("${prefix}_amrfinderplus.tsv"), emit: report
    tuple val("${task.process}"), val("amrfinderplus"), eval("amrfinder --version 2>&1 | sed -n '1p' || echo 'unknown'"), topic: versions, emit: versions_amrfinderplus


    script:
    def args = task.ext.args ?: "--plus"
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    gzip -cdf ${faa} > query.faa
    gzip -cdf ${fna} > query.fna
    gzip -cdf ${gff} > query.gff3

    if [ ! -s query.faa ] && [ ! -s query.fna ]; then
        printf 'Protein id\\tContig id\\tStart\\tStop\\tStrand\\tElement symbol\\tElement name\\tScope\\tType\\tSubtype\\tClass\\tSubclass\\tMethod\\t%% Coverage of reference\\t%% Identity to reference\\tClosest reference accession\\n' > ${prefix}_amrfinderplus.tsv
    else
        DB_DIR=\$(find -L ${db} -maxdepth 2 -type f -name "version.txt" -exec dirname {} \\; | sed -n '1p')
        DB_DIR=\${DB_DIR:-${db}}

        NEED_PROT_IDX=0
        NEED_DNA_IDX=0
        NEED_HMM_IDX=0
        if [ -f "\${DB_DIR}/AMRProt.fa" ] && [ ! -f "\${DB_DIR}/AMRProt.fa.pdb" ] && [ ! -f "\${DB_DIR}/AMRProt.fa.pin" ] && [ ! -f "\${DB_DIR}/AMRProt.pdb" ] && [ ! -f "\${DB_DIR}/AMRProt.pin" ]; then
            NEED_PROT_IDX=1
        fi
        for dna_fa in "\${DB_DIR}"/AMR_DNA-*.fa; do
            if [ -f "\${dna_fa}" ] && [ ! -f "\${dna_fa}.ndb" ] && [ ! -f "\${dna_fa}.nin" ]; then
                NEED_DNA_IDX=1
                break
            fi
        done
        if [ -f "\${DB_DIR}/AMR.LIB" ] && [ ! -f "\${DB_DIR}/AMR.LIB.h3m" ]; then
            NEED_HMM_IDX=1
        fi
        if [ "\${NEED_PROT_IDX}" -eq 1 ] || [ "\${NEED_DNA_IDX}" -eq 1 ] || [ "\${NEED_HMM_IDX}" -eq 1 ]; then
            cp -rL "\${DB_DIR}" amrfinder_db_work
            DB_DIR="amrfinder_db_work"
            if [ "\${NEED_PROT_IDX}" -eq 1 ]; then
                makeblastdb -in "\${DB_DIR}/AMRProt.fa" -dbtype prot -out "\${DB_DIR}/AMRProt.fa" >/dev/null
                makeblastdb -in "\${DB_DIR}/AMRProt.fa" -dbtype prot -out "\${DB_DIR}/AMRProt" >/dev/null
            fi
            if [ "\${NEED_DNA_IDX}" -eq 1 ]; then
                for dna_fa in "\${DB_DIR}"/AMR_DNA-*.fa; do
                    if [ -f "\${dna_fa}" ] && [ ! -f "\${dna_fa}.ndb" ] && [ ! -f "\${dna_fa}.nin" ]; then
                        makeblastdb -in "\${dna_fa}" -dbtype nucl -out "\${dna_fa}" >/dev/null
                    fi
                done
            fi
            if [ "\${NEED_HMM_IDX}" -eq 1 ]; then
                hmmpress -f "\${DB_DIR}/AMR.LIB" >/dev/null
            fi
        fi

        HAS_FAA=0
        HAS_FNA=0
        HAS_GFF=0
        if [ -s query.faa ]; then
            HAS_FAA=1
        fi
        if [ -s query.fna ]; then
            HAS_FNA=1
        fi
        if [ -s query.gff3 ] && grep -q '^[^#]' query.gff3; then
            HAS_GFF=1
        fi

        if [ "\${HAS_FAA}" -eq 1 ] && [ "\${HAS_FNA}" -eq 1 ] && [ "\${HAS_GFF}" -eq 1 ]; then
            INPUT_ARGS="--protein query.faa --nucleotide query.fna --gff query.gff3 --annotation_format standard"
        elif [ "\${HAS_FAA}" -eq 1 ]; then
            INPUT_ARGS="--protein query.faa"
        else
            INPUT_ARGS="--nucleotide query.fna"
        fi

        amrfinder \\
            ${args} \\
            \${INPUT_ARGS} \\
            --database "\${DB_DIR}" \\
            --threads ${task.cpus} \\
            --output ${prefix}_amrfinderplus.tsv
    fi
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    printf 'Protein id\\tContig id\\tStart\\tStop\\tStrand\\tElement symbol\\tElement name\\tScope\\tType\\tSubtype\\tClass\\tSubclass\\tMethod\\t%% Coverage of reference\\t%% Identity to reference\\tClosest reference accession\\n' > ${prefix}_amrfinderplus.tsv
    """
}
