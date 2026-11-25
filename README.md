# advena

**advena** (v2.0.0) is a modular Nextflow DSL2 pipeline for **mobileome discovery, cargo profiling, and comparative
genomics** of Integrative and Conjugative Elements (ICEs) and Integrative and Mobilizable Elements (IMEs) in _Bacillota_
genomes.

Built around a Nextflow DSL2 reimplementation of [ICEscreen](https://github.com/ICEscreen/ICEscreen) v1.3.3, **advena**
keeps core ICE/IME detection and element sequence extraction enabled by default while providing opt-in subworkflows
(`--run_*` flags or `--run_all`) for raw FASTA annotation, assembly/phylum QC, dereplicated AMR/stress profiling, eggNOG
Orthologous Group (OG) classification, anti-phage defense & BGC discovery, plasmid/prophage/IS/integron/`att` context,
and cohort clustering, MultiQC dashboards, and synteny plots.

## Features

- **Flexible input & pre-detection QC**: Accept annotated GenBank (`.gb`/`.gbk`/`.gbff[.gz]`), raw FASTA assemblies
  (auto-annotated with **Bakta**), or both; optional **QUAST** + **CheckM2** + **GTDB-Tk** QC/phylum gate (`--run_qc`,
  `--run_taxonomy_check`)
- **Core ICEscreen v1.3.3 + always-on element slicing**: Detects ICE/IME structures and slices each element into BED6,
  `.fna.gz`, `.faa.gz`, `.gff3.gz`, CDS manifest, and per-element `.gbk.gz` (`--me_flank_bp`)
- **Mobile cargo profiling**:
  - **AMRFinderPlus** (`--plus`) + **CARD RGI** (Strict + Perfect) merged via native locus-tag and >=80% reciprocal
    coordinate overlap dereplication (`DEREP_AMR_CARD`) preserving AMR, metal, biocide, heat/acid stress, and virulence
    genes alongside CARD ARO ontology terms (`--run_amr`)
  - **eggNOG-mapper** (DIAMOND mode) extracting hierarchical `Bacteria (@2)` and `Bacillota (@1239)` OGs, COG
    categories, KEGG, and PFAMs (`--run_eggnog`)
  - **DefenseFinder** (MacSyFinder + CasFinder models) for anti-phage defense and TA systems (`--run_defense`)
  - **antiSMASH-lite** secondary metabolite/bacteriocin BGC intersection (`--run_bgc`)
- **Complementary MGE context**: **MOB-suite** recon, **geNomad**, **ISEScan**, and **IntegronFinder**
  (`--run_mge_context`), plus `attL`/`attR` direct repeat and tRNA/housekeeping hotspot detection (`--run_att_sites`)
- **Cohort comparative genomics & reporting**: **skani** (`--small-genomes`) + **MMseqs2** element clustering
  (`--run_cohort_cluster`), pan-ME eggNOG OG presence/absence matrix, **MultiQC** dashboard (`--run_multiqc`), and
  **pyGenomeViz** SVG/HTML synteny alignments (`--run_synteny`)
- **Unified database caching & 4-tier local verification**: `--db_cache_dir` (`storeDir`), per-tool `--*_db` overrides,
  > 30 GB DB safety guards, `pytest` + `Ruff`, `nf-test`, and local GitHub Actions testing via `gh act`

## Requirements

- [Nextflow](https://www.nextflow.io/) >= 26.04.0 (Java >= 17)
- One of: Docker, Singularity, Apptainer, Podman, or Conda

## Quick Start

```bash
# Default core run (ICEscreen + element sequence extraction)
nextflow run andreassag/advena \
    --input assets/samplesheet.csv \
    --outdir results \
    -profile docker

# Enable cargo profiling + MGE context + cohort reporting
nextflow run andreassag/advena \
    --input assets/samplesheet.csv \
    --outdir results \
    --db_cache_dir /data/advena_dbs \
    --eggnog_db /data/eggnog_db \
    --run_amr --run_eggnog --run_defense --run_att_sites --run_cohort_cluster --run_multiqc --run_synteny \
    -profile docker
```

Verify with the bundled test datasets:

```bash
nextflow run andreassag/advena -profile test,docker
nextflow run andreassag/advena -profile test_full,docker
```

See the [documentation](https://andreassag.github.io/advena) for full parameter, database, and output references.

## Input

A CSV samplesheet with one row per genome (`genbank`, `fasta`, or both):

| Column    | Description                                                      |
| --------- | ---------------------------------------------------------------- |
| `sample`  | Unique sample identifier (no spaces)                             |
| `genbank` | Optional path to GenBank file (`.gb`, `.gbk`, `.gbff`, `.gz`)    |
| `fasta`   | Optional path to FASTA assembly (`.fa`, `.fna`, `.fasta`, `.gz`) |

```csv
sample,genbank,fasta
genome_A,/path/to/genome_A.gbff,
genome_B,,/path/to/genome_B.fna.gz
```

## Output (v2.0.0 Layout)

```text
<outdir>/
├── bakta/<sample>/                 # Bakta annotation (FASTA-only inputs)
├── qc/                             # QUAST, CheckM2, GTDB-Tk & qc_gate_summary.tsv
├── icescreen/<sample>/             # Core ICEscreen results + me_sequences/ slices
├── cargo/<sample>/                 # AMRFinderPlus, CARD RGI, DEREP, eggNOG, DefenseFinder, BGC, advena.gff.gz
├── mge_context/<sample>/           # MOB-suite, geNomad, ISEScan, IntegronFinder, att_sites, contextualized ME TSV
├── cohort/                         # skani, MMseqs2, me_clusters.tsv, og_presence_absence.tsv, synteny/
├── multiqc/multiqc_report.html     # Interactive cohort MultiQC dashboard
└── pipeline_info/                  # Execution reports & software_versions.yml
```

## Local Development & Pre-Push Verification

```bash
make lint        # Nextflow lint, Prettier, Ruff, JSON schemas, MkDocs strict
make test-unit   # pytest unit tests on bin/ Python helpers
make test-stub   # Default + --run_all DAG stubs and nf-test --tag stub
make test        # Containerized nf-test suite (--profile docker)
make test-ci     # Run GitHub Actions PR workflows locally via gh act
make preflight   # All-in-one pre-push gate
```

## Citations

See [CITATIONS.md](CITATIONS.md) for the complete list of tools and references used by **advena**.

## License

This project is licensed under the [GNU Affero General Public License v3.0](LICENSE).
