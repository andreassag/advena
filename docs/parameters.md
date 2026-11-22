# Parameters

## Input / output

| Parameter  | Type   | Default    | Description                 |
| ---------- | ------ | ---------- | --------------------------- |
| `--input`  | string | (required) | Path to the samplesheet CSV |
| `--outdir` | string | `results`  | Directory for output files  |

---

## ICEscreen

| Parameter  | Type   | Default                                                                                      | Description                                                                        |
| ---------- | ------ | -------------------------------------------------------------------------------------------- | ---------------------------------------------------------------------------------- |
| `--phylum` | string | `bacillota`                                                                                  | Taxonomic phylum; selects `mode/<phylum>.yml` (currently only `bacillota`)         |
| `--db_url` | string | `https://github.com/andreassag/icescreen/releases/download/v1.3.3/icescreen_database.tar.gz` | URL or local path to the ICEscreen database archive (`.tar.gz`)                    |
| `--db_dir` | string | `null`                                                                                       | Local directory containing a pre-extracted database or cache target for extraction |

---

## BLASTP

| Parameter                  | Type    | Default | Description                                  |
| -------------------------- | ------- | ------- | -------------------------------------------- |
| `--blastp_evalue`          | float   | `0.001` | E-value threshold for BLASTP hits            |
| `--blastp_max_target_seqs` | integer | `10`    | Maximum number of target sequences per query |

---

## Generic options

| Parameter               | Type           | Default               | Description                                                 |
| ----------------------- | -------------- | --------------------- | ----------------------------------------------------------- |
| `--publish_dir_mode`    | string         | `copy`                | Nextflow `publishDir` mode: `copy`, `symlink`, `move`, etc. |
| `--monochrome_logs`     | boolean        | `false`               | Disable ANSI colour in log output                           |
| `--help`                | boolean/string | `false`               | Display the help message (or parameter help) and exit       |
| `--version`             | boolean        | `false`               | Display the pipeline version and exit                       |
| `--trace_report_suffix` | string         | `yyyy-MM-dd_HH-mm-ss` | Timestamp suffix appended to `pipeline_info/` reports       |
