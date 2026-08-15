# Usage

## Samplesheet

advena takes a CSV samplesheet as input, with one genome per row:

```csv
sample,genbank
GCF_000009045,/data/GCF_000009045.gbff
GCF_003176835,/data/GCF_003176835.gbff
```

| Column    | Description                                           |
| --------- | ----------------------------------------------------- |
| `sample`  | Unique sample identifier used in all output filenames |
| `genbank` | Path to a GenBank file (`.gb` or `.gbff`)             |

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

# Conda (local)
-profile conda

# Test dataset with Docker
-profile test,docker
```

---

## Specifying the ICEscreen installation

advena reads the ICEscreen Python scripts and databases at runtime. There are three ways to provide the paths:

1. **Environment variable (recommended)**: Set `ICESCREEN_ROOT` before running the pipeline.

```bash
export ICESCREEN_ROOT=/opt/icescreen
nextflow run andreassag/advena --input samplesheet.csv --outdir results -profile docker
```

2. **Pipeline parameter**: Pass `--icescreen_root` directly.

```bash
nextflow run andreassag/advena \
  --input samplesheet.csv \
  --outdir results \
  --icescreen_root /opt/icescreen \
  -profile docker
```

3. **Override database path**: Use `--icescreen_db` to override database path.

```bash
nextflow run andreassag/advena \
  --input samplesheet.csv \
  --outdir results \
  --icescreen_db /data/icescreen_databases \
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

---

## Adjusting ME detection parameters

Mobile element segmentation thresholds can be adjusted for non-standard genomes:

```bash
nextflow run andreassag/advena \
    --input samplesheet.csv \
    --outdir results \
    --min_cds_between_segments 50 \
    --max_cds_for_ime_size 20 \
    -profile docker
```
