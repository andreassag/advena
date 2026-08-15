# advena

**advena** is a Nextflow pipeline for ICEscreen: detection and annotation of Integrative and Conjugative Elements (ICEs)
and Integrative and Mobilizable Elements (IMEs) in Bacillota genomes.

This pipeline is a Nextflow DSL2 reimplementation of the original [ICEscreen](https://github.com/ICEscreen/ICEscreen)
Snakemake pipeline, orchestrating the ICEscreen Python scripts and databases via Docker or Conda.

---

## Overview

```mermaid
graph LR
    A[GenBank input<br/>samplesheet.csv] --> B[GB_TO_FAA]
    B --> C[BLASTP_SEARCH<br/>×7 databases]
    B --> D[HMMSCAN_SP<br/>~20 profiles]
    C --> E[PROCESS_BLASTP]
    D --> F[PROCESS_HMMSCAN]
    E --> G[MERGE_SP]
    F --> G
    G --> H[FP_SCREENING]
    H --> I[REANNOT_SP]
    I --> J[DETECT_ME]
    I --> K[CREATE_ANNOTATIONS]
    J --> K
    J --> L[PACKAGE_RESULTS]
    K --> M[GFF3 / EMBL / GenBank]
```

---

## Features

- ICE/IME structure detection and boundary assembly using the original
  [ICEscreen](https://github.com/ICEscreen/ICEscreen) logic
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
