# Usage

## Samplesheet

advena takes a CSV samplesheet as input, with one genome per row:

```csv
sample,genbank
GCF_000009045.1,/data/GCF_000009045.1.gbff
GCF_003176835.1,/data/GCF_003176835.1.gbff
```

| Column    | Description                                                   |
| --------- | ------------------------------------------------------------- |
| `sample`  | Unique sample identifier (no spaces) used in output filenames |
| `genbank` | Path to a GenBank file (`.gb`, `.gbk`, or `.gbff`)            |

An example samplesheet is provided in
[`assets/samplesheet.csv`](https://github.com/andreassag/advena/blob/main/assets/samplesheet.csv).

---

## Basic usage

```bash
nextflow run andreassag/advena \
    --input samplesheet.csv \
    --outdir results \
    -profile docker
```

---

## Container engine profiles

Combine a container engine profile with an optional execution profile:

```bash
# Docker (local)
-profile docker

# Singularity on an HPC cluster
-profile singularity,hpc

# SLURM cluster with Singularity
-profile singularity,slurm

# Apptainer (local)
-profile apptainer

# Podman (local)
-profile podman

# Conda (local)
-profile conda

# Test dataset with Docker
-profile test,docker
```

---

## Database caching and custom archives

By default, advena downloads `icescreen_database.tar.gz` from `--db_url`, extracts it, builds protein BLAST indices, and
presses HMM profiles automatically.

1. **Cache or reuse an extracted database directory (`--db_dir`)**:

```bash
nextflow run andreassag/advena \
    --input samplesheet.csv \
    --outdir results \
    --db_dir /data/icescreen_db_cache \
    -profile docker
```

If `/data/icescreen_db_cache/icescreen_database` (or `/data/icescreen_db_cache` directly) already contains
`icescreen.conf`, extraction is skipped and the cached files are reused directly. Otherwise, `UNTAR` extracts the
archive into `--db_dir` for future runs.

2. **Use a custom archive URL or local `.tar.gz` file (`--db_url`)**:

```bash
nextflow run andreassag/advena \
    --input samplesheet.csv \
    --outdir results \
    --db_url /data/icescreen_database.tar.gz \
    -profile docker
```

---

## Resuming a run

Nextflow caches completed tasks. If a run is interrupted, resume from where it left off:

```bash
nextflow run andreassag/advena \
    --input samplesheet.csv \
    --outdir results \
    -profile docker \
    -resume
```

---

## Adjusting BLAST parameters

The BLASTP search sensitivity can be tuned with `--blastp_evalue` and `--blastp_max_target_seqs`:

```bash
nextflow run andreassag/advena \
    --input samplesheet.csv \
    --outdir results \
    --blastp_evalue 1e-5 \
    --blastp_max_target_seqs 5 \
    -profile docker
```
