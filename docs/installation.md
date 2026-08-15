# Installation

## Requirements

| Dependency       | Minimum version | Notes                                                   |
| ---------------- | --------------- | ------------------------------------------------------- |
| Nextflow         | 23.04.0         | `NXF_VER` env var can pin the version                   |
| Java             | 11              | Required by Nextflow                                    |
| Container engine | —               | One of Docker, Singularity, Apptainer, Podman, or Conda |
| ICEscreen        | —               | Must be installed and accessible via `ICESCREEN_ROOT`   |

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

advena requires access to the ICEscreen software environment. The recommended approach is Docker or Singularity with the
pre-built container.

=== "Docker"

    ```bash
    # Install Docker: https://docs.docker.com/engine/install/
    docker pull andreassag/icescreen-nf:latest
    ```

=== "Singularity / Apptainer"

    ```bash
    singularity pull icescreen-nf.sif docker://andreassag/icescreen-nf:latest
    ```

=== "Conda"

    ```bash
    # Install Mamba for faster environment resolution
    conda install -n base -c conda-forge mamba
    ```

---

## 3. Install ICEscreen

advena calls ICEscreen Python scripts and databases at runtime. Point the pipeline to an ICEscreen installation by
setting the `ICESCREEN_ROOT` environment variable or the `--icescreen_root` parameter.

```bash
# Clone the ICEscreen repository
git clone https://github.com/ICEscreen/ICEscreen.git /opt/icescreen

# Export the root path
export ICESCREEN_ROOT=/opt/icescreen
```

If `ICESCREEN_ROOT` is not set, advena falls back to `~/repo/icescreen` then `/opt/icescreen`. If the databases are
stored separately, override the database path with `--icescreen_db`.

---

## 4. Run advena

```bash
nextflow run andreassag/advena \
    --input samplesheet.csv \
    --outdir results \
    -profile docker
```

To pin a specific version:

```bash
NXF_VER=23.10.0 nextflow run andreassag/advena -r 1.0.0 \
    --input samplesheet.csv \
    --outdir results \
    -profile docker
```
