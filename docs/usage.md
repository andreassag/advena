# Usage

## Samplesheet

advena accepts a CSV samplesheet with one genome per row. Provide `genbank`, `fasta`, or both:

```csv
sample,genbank,fasta
GCF_000009045.1,/data/GCF_000009045.1.gbff,
isolate_raw,,/data/isolate_raw.fna.gz
```

| Column    | Description                                                       |
| --------- | ----------------------------------------------------------------- |
| `sample`  | Unique sample identifier (no spaces) used in output filenames     |
| `genbank` | Optional GenBank file (`.gb`, `.gbk`, or `.gbff`, optional `.gz`) |
| `fasta`   | Optional FASTA assembly (`.fa`, `.fna`, `.fasta`, optional `.gz`) |

When only `fasta` is supplied, advena automatically annotates the assembly using **Bakta** before running ICEscreen.

---

## Running Core vs. Extended Modules

By default, advena runs core ICEscreen detection and element sequence extraction (`EXTRACT_ME_SEQUENCES`). External
database stages are opt-in:

```bash
# Default core run
nextflow run andreassag/advena \
    --input samplesheet.csv \
    --outdir results \
    -profile docker

# Cargo profiling (AMRFinderPlus + CARD RGI + eggNOG + DefenseFinder) + cohort reporting
nextflow run andreassag/advena \
    --input samplesheet.csv \
    --outdir results \
    --db_cache_dir /data/advena_dbs \
    --eggnog_db /data/eggnog_db \
    --run_amr --run_eggnog --run_defense \
    --run_att_sites --run_cohort_cluster --run_multiqc --run_synteny \
    -profile docker

# Enable every optional module at once
nextflow run andreassag/advena \
    --input samplesheet.csv \
    --outdir results \
    --db_cache_dir /data/advena_dbs \
    --eggnog_db /data/eggnog_db \
    --gtdbtk_db /data/gtdbtk_db \
    --run_all \
    -profile docker
```

---

## Database Caching & Large DB Safety Guards

Pass `--db_cache_dir /path/to/cache` to automatically download and reuse small/medium reference databases via Nextflow
`storeDir` (ICEscreen, Bakta-light, CheckM2, AMRFinderPlus, CARD, DefenseFinder, MOB-suite, and geNomad).

Any database can be overridden with an explicit local path (`--amrfinderplus_db`, `--card_db`, `--eggnog_db`,
`--defensefinder_db`, `--bakta_db`, `--checkm2_db`, `--gtdbtk_db`, `--mobsuite_db`, `--genomad_db`, `--antismash_db`).

For databases exceeding 30 GB (**GTDB-Tk**, **eggNOG-mapper**, and **Bakta full**), advena requires an explicit `--*_db`
path unless `--auto_download_large_dbs true` is set.

---

## Container Engine & Execution Profiles

```bash
-profile docker               # Docker (local)
-profile singularity,hpc      # Singularity on an HPC cluster
-profile singularity,slurm    # SLURM cluster with Singularity
-profile apptainer            # Apptainer
-profile podman               # Podman
-profile conda                # Conda (local)
-profile test,docker          # Core Bacillota golden test dataset
-profile test_full,docker     # Offline integration test using miniature fixture DBs
```
