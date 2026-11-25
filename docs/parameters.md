# Parameters

## Input / output

| Parameter  | Type   | Default    | Description                 |
| ---------- | ------ | ---------- | --------------------------- |
| `--input`  | string | (required) | Path to the samplesheet CSV |
| `--outdir` | string | `results`  | Directory for output files  |

---

## ICEscreen & BLASTP

| Parameter                  | Type    | Default                                                                                      | Description                                              |
| -------------------------- | ------- | -------------------------------------------------------------------------------------------- | -------------------------------------------------------- |
| `--phylum`                 | string  | `bacillota`                                                                                  | Taxonomic phylum (`bacillota`)                           |
| `--db_url`                 | string  | `https://github.com/andreassag/icescreen/releases/download/v1.3.3/icescreen_database.tar.gz` | URL or local path to ICEscreen DB archive (`.tar.gz`)    |
| `--db_dir`                 | string  | `null`                                                                                       | Pre-extracted ICEscreen DB directory or extraction cache |
| `--blastp_evalue`          | float   | `0.001`                                                                                      | E-value threshold for BLASTP hits                        |
| `--blastp_max_target_seqs` | integer | `10`                                                                                         | Maximum number of target sequences per BLASTP query      |

---

## Optional module toggles

| Parameter              | Type    | Default | Description                                                                   |
| ---------------------- | ------- | ------- | ----------------------------------------------------------------------------- |
| `--run_all`            | boolean | `false` | Convenience preset enabling all optional QC, cargo, MGE, and cohort modules   |
| `--run_qc`             | boolean | `false` | Run QUAST and CheckM2 assembly quality and completeness checks                |
| `--run_taxonomy_check` | boolean | `false` | Run GTDB-Tk taxonomic classification for Bacillota phylum verification        |
| `--run_amr`            | boolean | `false` | Run AMRFinderPlus (`--plus`), CARD RGI, and native `DEREP_AMR_CARD`           |
| `--run_eggnog`         | boolean | `false` | Run eggNOG-mapper (DIAMOND) and extract `@2` Bacteria + `@1239` Bacillota OGs |
| `--run_defense`        | boolean | `false` | Run DefenseFinder (MacSyFinder + CasFinder models)                            |
| `--run_bgc`            | boolean | `false` | Run antiSMASH-lite on GenBank and intersect BGC regions with ICE/IMEs         |
| `--run_mge_context`    | boolean | `false` | Run MOB-suite recon, geNomad, ISEScan, and IntegronFinder                     |
| `--run_att_sites`      | boolean | `false` | Detect `attL`/`attR` direct repeats and tRNA/housekeeping targets             |
| `--run_cohort_cluster` | boolean | `false` | Run `skani` (`--small-genomes`) + `MMseqs2` cohort ICE/IME clustering         |
| `--run_multiqc`        | boolean | `false` | Generate interactive cohort MultiQC dashboard                                 |
| `--run_synteny`        | boolean | `false` | Render `pyGenomeViz` SVG and HTML cluster synteny alignments                  |

---

## Database options

| Parameter                   | Type    | Default | Description                                                      |
| --------------------------- | ------- | ------- | ---------------------------------------------------------------- |
| `--db_cache_dir`            | string  | `null`  | Unified `storeDir` cache for auto-downloaded reference databases |
| `--auto_download_large_dbs` | boolean | `false` | Allow auto-download of >30 GB DBs (GTDB-Tk, eggNOG, Bakta full)  |
| `--bakta_db`                | string  | `null`  | Pre-downloaded Bakta database directory                          |
| `--bakta_db_type`           | string  | `light` | Bakta auto-download flavour (`light` or `full`)                  |
| `--checkm2_db`              | string  | `null`  | Pre-downloaded CheckM2 DIAMOND database                          |
| `--gtdbtk_db`               | string  | `null`  | Pre-downloaded GTDB-Tk reference data directory                  |
| `--amrfinderplus_db`        | string  | `null`  | Pre-downloaded AMRFinderPlus database directory                  |
| `--card_db`                 | string  | `null`  | Pre-downloaded CARD database directory containing `card.json`    |
| `--eggnog_db`               | string  | `null`  | Pre-downloaded eggNOG-mapper database directory                  |
| `--defensefinder_db`        | string  | `null`  | Pre-downloaded DefenseFinder models directory                    |
| `--mobsuite_db`             | string  | `null`  | Pre-downloaded MOB-suite database directory                      |
| `--genomad_db`              | string  | `null`  | Pre-downloaded geNomad database directory                        |
| `--antismash_db`            | string  | `null`  | Pre-downloaded antiSMASH database directory                      |

---

## QC, cargo, MGE context, and cohort options

| Parameter                            | Type    | Default    | Description                                                                  |
| ------------------------------------ | ------- | ---------- | ---------------------------------------------------------------------------- |
| `--strict_qc`                        | boolean | `false`    | Exclude samples failing CheckM2 thresholds instead of warning only           |
| `--min_completeness`                 | float   | `90.0`     | Minimum CheckM2 completeness (%)                                             |
| `--max_contamination`                | float   | `5.0`      | Maximum CheckM2 contamination (%)                                            |
| `--filter_non_bacillota`             | boolean | `false`    | Exclude non-Bacillota samples instead of warning only                        |
| `--me_flank_bp`                      | integer | `0`        | Flanking bp upstream/downstream when slicing ICE/IME intervals               |
| `--cargo_scope`                      | string  | `elements` | Cargo profiling scope: `elements` (extracted MEs only) or `genome`           |
| `--amr_derep_overlap`                | float   | `0.8`      | Minimum reciprocal coordinate overlap for AMRFinderPlus + CARD dereplication |
| `--rgi_include_loose`                | boolean | `false`    | Include CARD RGI `Loose` hits in addition to `Perfect` and `Strict`          |
| `--eggnog_tax_scope`                 | string  | `1239`     | eggNOG-mapper `--tax_scope` (default `1239` Bacillota/Firmicutes)            |
| `--att_flank_bp`                     | integer | `5000`     | Flanking window (bp) searched around boundary integrases for `attL`/`attR`   |
| `--me_ani_threshold`                 | float   | `95.0`     | Minimum skani ANI (%) to cluster homologous ICEs/IMEs                        |
| `--me_af_threshold`                  | float   | `70.0`     | Minimum skani aligned fraction (%) to cluster homologous ICEs/IMEs           |
| `--synteny_max_elements_per_cluster` | integer | `15`       | Maximum ICE/IME slices plotted per cluster in `pyGenomeViz` figures          |

---

## Generic options

| Parameter               | Type           | Default               | Description                                                 |
| ----------------------- | -------------- | --------------------- | ----------------------------------------------------------- |
| `--publish_dir_mode`    | string         | `copy`                | Nextflow `publishDir` mode: `copy`, `symlink`, `move`, etc. |
| `--monochrome_logs`     | boolean        | `false`               | Disable ANSI colour in log output                           |
| `--help`                | boolean/string | `false`               | Display the help message (or parameter help) and exit       |
| `--version`             | boolean        | `false`               | Display the pipeline version and exit                       |
| `--trace_report_suffix` | string         | `yyyy-MM-dd_HH-mm-ss` | Timestamp suffix appended to `pipeline_info/` reports       |
