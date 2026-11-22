# advena

**advena** is a Nextflow pipeline for ICEscreen: detection and annotation of Integrative and Conjugative Elements (ICEs)
and Integrative and Mobilizable Elements (IMEs) in Bacillota genomes.

This pipeline is a Nextflow DSL2 reimplementation of the original [ICEscreen](https://github.com/ICEscreen/ICEscreen)
Snakemake pipeline, combining standard [`nf-core/modules`](https://github.com/nf-core/modules) with containerized
[ICEscreen](https://github.com/andreassag/icescreen) Python tools (`ghcr.io/andreassag/icescreen:v1.3.3`).

---

## Overview

```mermaid
graph LR
    A[GenBank input<br/>samplesheet.csv] --> B[GB_TO_FAA]
    DB[db_url / db_dir] --> P[PREPARE_DATABASE<br/>UNTAR + MAKEBLASTDB + HMMPRESS]
    P --> B
    B --> C[BLAST_BLASTP<br/>x7 databases]
    P --> C
    B --> D[HMMSCAN_SP<br/>22 profiles]
    P --> D
    C --> E[PROCESS_BLASTP]
    D --> F[PROCESS_HMMSCAN]
    E --> G[MERGE_SP]
    F --> G
    G --> H1[EXTRACT_SP_FAA]
    B --> H1
    H1 --> H2[HMMSCAN_FP<br/>3 FP profiles]
    P --> H2
    G --> H3[PROCESS_HMMSCAN_FP]
    H2 --> H3
    H3 --> I[REANNOT_SP]
    E --> I
    I --> J[DETECT_ME]
    A --> J
    I --> K[CREATE_ANNOTATIONS]
    J --> K
    J --> L[PACKAGE_RESULTS]
    K --> M[GFF3 / EMBL / GenBank]
```

---

## Features

- ICE/IME structure detection and boundary assembly using the [ICEscreen](https://github.com/andreassag/icescreen)
  v1.3.3 logic
- Automatic database retrieval and indexing via `nf-core/untar`, `nf-core/blast/makeblastdb`, and
  `nf-core/hmmer/hmmpress`
- Parameter and samplesheet validation via `nf-schema`
- Support for Docker, Singularity, Apptainer, Podman, and Conda

---

## Quick start

```bash
nextflow run andreassag/advena \
    --input samplesheet.csv \
    --outdir results \
    -profile docker
```

See [Installation](installation.md) and [Usage](usage.md) for full details.
