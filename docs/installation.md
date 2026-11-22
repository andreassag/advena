# Installation

## Requirements

| Dependency       | Minimum version | Notes                                                   |
| ---------------- | --------------- | ------------------------------------------------------- |
| Nextflow         | 25.10.0         | `NXF_VER` env var can pin the version                   |
| Java             | 17              | Required by Nextflow and `nf-schema`                    |
| Container engine | —               | One of Docker, Singularity, Apptainer, Podman, or Conda |

---

## 1. Install Nextflow

```bash
curl -s https://get.nextflow.io | bash
# or via conda
conda install -c bioconda nextflow
```

Add the `nextflow` binary to your `PATH`, or move it to a location already on it.

---

## 2. Install a container engine

advena uses containerized software environments for ICEscreen, BLAST+, HMMER, and GNU coreutils. No separate local
ICEscreen installation is required; databases are automatically fetched from `--db_url` (or reused from `--db_dir`).

=== "Docker"

    ```bash
    # Install Docker: https://docs.docker.com/engine/install/
    docker pull ghcr.io/andreassag/icescreen:v1.3.3
    ```

=== "Singularity / Apptainer"

    ```bash
    singularity pull icescreen-v1.3.3.sif docker://ghcr.io/andreassag/icescreen:v1.3.3
    ```

=== "Conda"

    ```bash
    # Install Mamba for faster environment resolution
    conda install -n base -c conda-forge mamba
    ```

---

## 3. Run advena

Verify your installation using the bundled test dataset:

```bash
nextflow run andreassag/advena \
    -profile test,docker
```

Run on your own data:

```bash
nextflow run andreassag/advena \
    --input samplesheet.csv \
    --outdir results \
    -profile docker
```

To pin a specific pipeline version:

```bash
NXF_VER=25.10.0 nextflow run andreassag/advena -r 1.0.0 \
    --input samplesheet.csv \
    --outdir results \
    -profile docker
```
