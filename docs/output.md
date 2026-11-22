# Output

Results are written to `<outdir>/ICEscreen_results/<sample>/` for each sample, with execution metadata in
`<outdir>/pipeline_info/`:

```text
<outdir>/
├── ICEscreen_results/
│   └── <sample>/
│       ├── detected_mobile_elements/
│       │   ├── <sample>_detected_ME.tsv
│       │   ├── <sample>_detected_ME.summary
│       │   ├── <sample>_detected_SP_withMEIds.tsv
│       │   └── standard_genome_annotation_formats/
│       │       ├── <sample>_icescreen.gff.gz
│       │       ├── <sample>_icescreen.embl.gz
│       │       ├── <sample>_icescreen.gb.gz
│       │       ├── <sample>_source.fa.gz
│       │       └── <sample>_source.gff.gz
│       ├── detected_signature_proteins/
│       │   └── <sample>_detected_SP.tsv
│       ├── param.conf.gz
│       └── tmp_intermediate_files.tar.gz
└── pipeline_info/
    ├── execution_report_<timestamp>.html
    ├── execution_timeline_<timestamp>.html
    ├── execution_trace_<timestamp>.txt
    ├── pipeline_dag_<timestamp>.html
    └── software_versions.yml
```

---

## Detected mobile elements

### `<sample>_detected_ME.tsv`

Tab-separated table of detected ICE and IME structures. Each row represents one mobile element, with columns for genomic
coordinates, element type, boundary signature proteins, and associated CDS.

### `<sample>_detected_ME.summary`

Human-readable summary of detection results: counts of detected ICEs, IMEs, and unclassified elements per replicon.

### `<sample>_detected_SP_withMEIds.tsv`

Signature protein table annotated with the mobile element ID to which each protein was assigned.

---

## Genome annotation formats

All annotation files are gzip-compressed.

| File                         | Format  | Contents                                     |
| ---------------------------- | ------- | -------------------------------------------- |
| `<sample>_icescreen.gff.gz`  | GFF3    | ICE/IME features overlaid on the source GFF3 |
| `<sample>_icescreen.embl.gz` | EMBL    | ICE/IME features in EMBL flat-file format    |
| `<sample>_icescreen.gb.gz`   | GenBank | ICE/IME features in GenBank flat-file format |
| `<sample>_source.fa.gz`      | FASTA   | Genome nucleotide sequence                   |
| `<sample>_source.gff.gz`     | GFF3    | Source genome annotation (CDS features)      |

---

## Detected signature proteins

### `<sample>_detected_SP.tsv`

Final table of signature proteins after false positive removal and XerS re-annotation. Each row is one detected
signature protein with its type, locus tag, replicon, and position.

---

## Run archive

### `tmp_intermediate_files.tar.gz`

Archive of all intermediate analysis files, organised to mirror the original ICEscreen directory structure under
`results/<sample>/`. Includes BLASTP output per database, HMM scan output per profile, filtered results, merged SP
tables, FP screening results, and re-annotation files.

### `param.conf.gz`

Gzip-compressed parameter summary combining the ICEscreen detection log (search thresholds, counts of detected elements)
with the full detection mode configuration (from the ICEscreen mode YAML file). Mirrors the `param.conf` output of the
original Snakemake pipeline.

---

## Pipeline info

Nextflow execution reports and collated tool versions written to `pipeline_info/`:

| File                           | Contents                                |
| ------------------------------ | --------------------------------------- |
| `execution_report_<ts>.html`   | Per-task resource usage and status      |
| `execution_timeline_<ts>.html` | Gantt chart of task execution timing    |
| `execution_trace_<ts>.txt`     | Tab-separated task metadata             |
| `pipeline_dag_<ts>.html`       | Directed acyclic graph of the workflow  |
| `software_versions.yml`        | Versions of all tools used across steps |
