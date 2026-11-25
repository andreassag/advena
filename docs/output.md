# Output (v2.0.0)

Results follow the `nf-core` per-stage top-level directory layout:

```text
<outdir>/
├── bakta/<sample>/                              # Upstream assembly annotation (FASTA-only inputs)
├── qc/
│   ├── checkm2/<sample>/                        # Completeness & contamination (--run_qc)
│   ├── quast/<sample>/                          # Assembly contiguity metrics (--run_qc)
│   ├── gtdbtk/<sample>/                         # Taxonomic classification (--run_taxonomy_check)
│   └── qc_gate_summary.tsv                      # Unified QC & Bacillota phylum gate table
├── icescreen/<sample>/
│   ├── detected_mobile_elements/                # *_detected_ME.tsv, .summary, *_detected_SP_withMEIds.tsv
│   │   └── standard_genome_annotation_formats/  # *_icescreen.{gff,embl,gb}.gz, *_source.{fa,gff}.gz
│   ├── detected_signature_proteins/             # *_detected_SP.tsv
│   ├── me_sequences/                            # *_ME_regions.{bed,fna.gz}, *_ME_cds.{faa.gz,gff3.gz,manifest.tsv}, slices/*.gbk.gz
│   ├── param.conf.gz
│   └── tmp_intermediate_files.tar.gz
├── cargo/<sample>/
│   ├── amrfinderplus/                           # Raw AMRFinderPlus (--plus) TSV
│   ├── card_rgi/                                # Raw CARD RGI main TXT & JSON
│   ├── amr_dereplicated/                        # <sample>_amr_card_derep.tsv (merged non-redundant calls)
│   ├── eggnog/                                  # .emapper.annotations + <sample>_cargo_ogs.tsv
│   ├── defensefinder/                           # DefenseFinder systems & genes TSVs
│   ├── antismash/                               # antiSMASH-lite outputs + <sample>_bgc_me_overlaps.tsv
│   ├── <sample>_me_cargo_profile.tsv            # Master joined locus-level cargo annotation table
│   └── <sample>_advena.gff.gz                   # Enriched GFF3 with ICE/IME + AMR/CARD + OG + Defense + BGC tags
├── mge_context/<sample>/
│   ├── mobsuite/                                # MOB-recon contig reports & plasmid reconstructions
│   ├── genomad/                                 # geNomad plasmid & virus/prophage summaries
│   ├── isescan/                                 # Insertion sequence intervals & families
│   ├── integronfinder/                          # Integron & attC cassette calls
│   ├── att_sites/                               # <sample>_att_sites.tsv (attL/attR direct repeats & tRNA targets)
│   └── <sample>_detected_ME_contextualized.tsv  # ICE/IME calls annotated with plasmid/prophage/IS/att flags
├── cohort/
│   ├── skani/                                   # Pairwise element ANI matrix
│   ├── mmseqs/                                  # Pan-ME protein clusters
│   ├── me_clusters.tsv                          # Cross-genome ICE/IME cluster assignments
│   ├── og_presence_absence.tsv                  # Cohort-wide eggNOG OG matrix (CORE_BACKBONE / SHELL_CARGO / CLOUD_ACCESSORY)
│   ├── cohort_me_cargo_summary.tsv              # Cohort-wide element + cargo summary table
│   └── synteny/                                 # pyGenomeViz SVG + HTML locus alignments per cluster
├── multiqc/
│   └── multiqc_report.html                      # Interactive cohort dashboard (--run_multiqc)
└── pipeline_info/
    ├── execution_report_<timestamp>.html
    ├── execution_timeline_<timestamp>.html
    ├── execution_trace_<timestamp>.txt
    ├── pipeline_dag_<timestamp>.html
    └── software_versions.yml
```

---

## Core ICEscreen & Extracted Mobile Elements (`icescreen/<sample>/`)

- `detected_mobile_elements/`: `*_detected_ME.tsv`, `*_detected_ME.summary`, `*_detected_SP_withMEIds.tsv`, plus
  gzip-compressed GFF3, EMBL, and GenBank annotations
- `detected_signature_proteins/`: `*_detected_SP.tsv`
- `me_sequences/`: Always-on interval slices (`*_ME_regions.bed`, `*_ME_regions.fna.gz`, `*_ME_cds.faa.gz`,
  `*_ME_cds.gff3.gz`, `*_ME_cds_manifest.tsv`, and `slices/<sample>__<ME_ID>.gbk.gz`)

---

## Mobile Cargo Profiling (`cargo/<sample>/`)

- `amr_dereplicated/<sample>_amr_card_derep.tsv`: Non-redundant AMR, stress, metal, biocide, and virulence calls merging
  NCBI AMRFinderPlus symbols/subtypes with CARD ARO accessions, drug classes, and resistance mechanisms
- `eggnog/<sample>_cargo_ogs.tsv`: Hierarchical `root_cog` (`@1`), `bacteria_og` (`@2`), `bacillota_og` (`@1239`),
  `narrowest_og`, `cog_category`, KEGG/PFAM annotations, and `cargo_role` (`ICE_BACKBONE_SP`, `CARGO_CHARACTERIZED_OG`,
  `CARGO_HYPOTHETICAL_OG`)
- `<sample>_me_cargo_profile.tsv` & `<sample>_advena.gff.gz`: Master locus-level table and enriched GFF3 overlaying all
  enabled cargo layers

---

## Complementary MGE Context (`mge_context/<sample>/`)

- `<sample>_detected_ME_contextualized.tsv`: Extends `*_detected_ME.tsv` with `plasmid_borne_warning`,
  `mobsuite_molecule_type`, `genomad_plasmid_score`, `composite_prophage_ids`, `flanking_is_families`,
  `internal_is_count`, `integron_ids`, `att_target_feature`, `att_repeat_length`, and `att_repeat_seq`

---

## Cohort Comparative Genomics & MultiQC (`cohort/`, `multiqc/`)

- `cohort/me_clusters.tsv`: Connected-component ICE/IME clusters (`CLUSTER_001`, ...) combining `skani` ANI/AF and
  `MMseqs2` protein-cluster Jaccard similarity
- `cohort/og_presence_absence.tsv`: Pan-ICE/IME eggNOG OG matrix distinguishing conserved backbone OGs from accessory
  cargo
- `cohort/synteny/<cluster>_synteny.{svg,html}`: Gene-track synteny alignments colored by functional layer
- `multiqc/multiqc_report.html`: Interactive cohort dashboard combining QUAST, CheckM2, Bakta, GTDB-Tk, and advena
  custom-content panels
