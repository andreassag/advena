#!/usr/bin/env python3
"""Build cohort ME cargo summary, pan-ME eggNOG OG presence/absence matrix, and MultiQC panels."""

from __future__ import annotations

import argparse
import csv
import sys
from collections import Counter
from pathlib import Path

from gbk_utils import open_text

COHORT_SUMMARY_COLUMNS = [
    "sample",
    "me_id",
    "contig",
    "category",
    "iceline_format",
    "ice_superfamily",
    "ice_family",
    "start",
    "end",
    "length_bp",
    "me_cluster_id",
    "qc_status",
    "is_bacillota",
    "plasmid_borne_warning",
    "composite_prophage_ids",
    "flanking_is_families",
    "att_target_feature",
    "n_cds",
    "n_amr_genes",
    "amr_genes",
    "n_defense_systems",
    "defense_systems",
    "bgc_types",
    "n_characterized_ogs",
    "n_hypothetical_ogs",
]


def _load_qc_by_sample(qc_path: Path | None) -> dict[str, dict[str, str]]:
    res: dict[str, dict[str, str]] = {}
    if qc_path is None or not qc_path.exists() or qc_path.stat().st_size == 0:
        return res
    with open_text(qc_path, "rt") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            s = (row.get("sample") or "").strip()
            if s:
                res[s] = row
    return res


def _load_clusters(cluster_path: Path | None) -> dict[tuple[str, str], str]:
    res: dict[tuple[str, str], str] = {}
    if cluster_path is None or not cluster_path.exists() or cluster_path.stat().st_size == 0:
        return res
    with open_text(cluster_path, "rt") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            s = (row.get("sample") or "").strip()
            me = (row.get("me_id") or "").strip()
            cid = (row.get("me_cluster_id") or "-").strip()
            if s and me:
                res[(s, me)] = cid
    return res


def _sample_from_me_filename(path: Path) -> str:
    name = path.name
    for suffix in ("_detected_ME_contextualized.tsv", "_detected_ME.tsv"):
        if name.endswith(suffix):
            return name[: -len(suffix)]
    return path.stem


def build_cohort_outputs(
    detected_me_paths: list[Path],
    contextualized_me_paths: list[Path],
    cargo_profile_paths: list[Path],
    qc_summary_path: Path | None,
    clusters_path: Path | None,
    out_dir: Path,
) -> tuple[int, int]:
    """Generate cohort summary TSV, OG presence/absence matrix, and MultiQC custom panels."""
    out_dir.mkdir(parents=True, exist_ok=True)
    qc_by_sample = _load_qc_by_sample(qc_summary_path)
    cluster_by_elem = _load_clusters(clusters_path)

    # Prefer contextualized ME tables when provided for a sample
    me_files_by_sample: dict[str, Path] = {}
    for p in detected_me_paths:
        if p.exists() and p.stat().st_size > 0:
            me_files_by_sample[_sample_from_me_filename(p)] = p
    for p in contextualized_me_paths:
        if p.exists() and p.stat().st_size > 0:
            me_files_by_sample[_sample_from_me_filename(p)] = p

    # Load cargo profile rows grouped by (sample, me_id) and build OG matrix structures
    cargo_by_elem: dict[tuple[str, str], list[dict[str, str]]] = {}
    og_counts: dict[str, Counter[str]] = {}
    og_meta: dict[str, dict[str, str]] = {}
    element_col_keys: list[str] = []
    cog_by_sample: dict[str, Counter[str]] = {}
    amr_by_sample: dict[str, Counter[str]] = {}

    for p in cargo_profile_paths:
        if not p.exists() or p.stat().st_size == 0:
            continue
        with open_text(p, "rt") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                sample = (row.get("sample") or "").strip()
                me_id = (row.get("me_id") or "").strip()
                if not sample or not me_id or me_id == "GENOME_BG":
                    continue
                cargo_by_elem.setdefault((sample, me_id), []).append(row)
                col_key = f"{sample}|{me_id}"
                if col_key not in element_col_keys:
                    element_col_keys.append(col_key)

                # Track COG distribution
                cog_cat = (row.get("cog_category") or "-").strip()
                if cog_cat and cog_cat != "-":
                    for ch in cog_cat:
                        if ch.isalpha():
                            cog_by_sample.setdefault(sample, Counter())[ch.upper()] += 1

                # Track dereplicated AMR/CARD genes
                amr_gene = (row.get("amr_gene") or "-").strip()
                if amr_gene and amr_gene != "-":
                    for g in amr_gene.split(";"):
                        if g and g != "-":
                            amr_by_sample.setdefault(sample, Counter())[g] += 1

                # Track hierarchical OGs (prefer bacillota_og @1239, fall back to bacteria_og @2)
                bac_og = (row.get("eggnog_bacillota_og") or "-").strip()
                bact_og = (row.get("eggnog_bacteria_og") or "-").strip()
                og_id = bac_og if bac_og != "-" else bact_og
                og_level = "Bacillota@1239" if bac_og != "-" else "Bacteria@2"
                if og_id and og_id != "-":
                    og_counts.setdefault(og_id, Counter())[col_key] += 1
                    if og_id not in og_meta:
                        og_meta[og_id] = {
                            "og_level": og_level,
                            "cog_category": cog_cat,
                            "description": (row.get("eggnog_description") or "-").strip(),
                        }

    summary_rows: list[dict[str, str]] = []
    sample_me_counts: dict[str, Counter[str]] = {}
    mge_context_counts: dict[str, Counter[str]] = {}

    for sample, me_path in sorted(me_files_by_sample.items()):
        sample_me_counts.setdefault(sample, Counter())
        mge_context_counts.setdefault(sample, Counter())
        qc = qc_by_sample.get(sample, {})
        with open_text(me_path, "rt") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                me_id = (row.get("ICE_IME_id") or "").strip()
                if not me_id or me_id == "-":
                    continue
                col_key = f"{sample}|{me_id}"
                if col_key not in element_col_keys:
                    element_col_keys.append(col_key)

                cat = (row.get("Category_of_element") or "Unknown").strip()
                try:
                    start = int(float(row.get("Start_of_most_upstream_SP") or 0))
                    end = int(float(row.get("Stop_of_most_downstream_SP") or 0))
                except ValueError:
                    start, end = 0, 0
                length_bp = max(0, abs(end - start) + 1) if (start > 0 and end > 0) else 0

                cat_low = cat.lower()
                if "complete ice" in cat_low:
                    sample_me_counts[sample]["Complete_ICE"] += 1
                elif "complete ime" in cat_low:
                    sample_me_counts[sample]["Complete_IME"] += 1
                else:
                    sample_me_counts[sample]["Partial_or_Other"] += 1

                plasmid_warn = (row.get("plasmid_borne_warning") or "NA").strip()
                prophage = (row.get("composite_prophage_ids") or "-").strip()
                flanking_is = (row.get("flanking_is_families") or "-").strip()
                att_target = (row.get("att_target_feature") or "-").strip()
                if plasmid_warn == "TRUE":
                    mge_context_counts[sample]["Plasmid_borne_MEs"] += 1
                if prophage != "-":
                    mge_context_counts[sample]["Composite_prophage_MEs"] += 1
                if flanking_is != "-":
                    mge_context_counts[sample]["IS_flanked_MEs"] += 1
                if att_target != "-":
                    mge_context_counts[sample]["att_site_resolved_MEs"] += 1

                cds_rows = cargo_by_elem.get((sample, me_id), [])
                amr_genes = sorted(
                    {
                        g
                        for r in cds_rows
                        for g in (r.get("amr_gene") or "-").split(";")
                        if g and g != "-"
                    }
                )
                def_systems = sorted(
                    {
                        d
                        for r in cds_rows
                        for d in (r.get("defense_system") or "-").split(";")
                        if d and d != "-"
                    }
                )
                bgcs = sorted(
                    {
                        b
                        for r in cds_rows
                        for b in (r.get("bgc_type") or "-").split(";")
                        if b and b != "-"
                    }
                )
                n_char = sum(1 for r in cds_rows if r.get("cargo_role") == "CARGO_CHARACTERIZED_OG")
                n_hyp = sum(1 for r in cds_rows if r.get("cargo_role") == "CARGO_HYPOTHETICAL_OG")

                summary_rows.append(
                    {
                        "sample": sample,
                        "me_id": me_id,
                        "contig": (row.get("Genome_accession") or "-").strip(),
                        "category": cat,
                        "iceline_format": (row.get("ICEline_format") or "-").strip(),
                        "ice_superfamily": (
                            row.get("ICE_consensus_superfamily_SP_conj_module") or "-"
                        ).strip(),
                        "ice_family": (
                            row.get("ICE_consensus_family_SP_conj_module") or "-"
                        ).strip(),
                        "start": str(start),
                        "end": str(end),
                        "length_bp": str(length_bp),
                        "me_cluster_id": cluster_by_elem.get((sample, me_id), "-"),
                        "qc_status": qc.get("qc_status", "NA"),
                        "is_bacillota": qc.get("is_bacillota", "NA"),
                        "plasmid_borne_warning": plasmid_warn,
                        "composite_prophage_ids": prophage,
                        "flanking_is_families": flanking_is,
                        "att_target_feature": att_target,
                        "n_cds": str(len(cds_rows)),
                        "n_amr_genes": str(len(amr_genes)),
                        "amr_genes": ";".join(amr_genes) or "-",
                        "n_defense_systems": str(len(def_systems)),
                        "defense_systems": ";".join(def_systems) or "-",
                        "bgc_types": ";".join(bgcs) or "-",
                        "n_characterized_ogs": str(n_char),
                        "n_hypothetical_ogs": str(n_hyp),
                    }
                )

    # 1. cohort_me_cargo_summary.tsv
    with open(out_dir / "cohort_me_cargo_summary.tsv", "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=COHORT_SUMMARY_COLUMNS, delimiter="\t")
        writer.writeheader()
        writer.writerows(summary_rows)

    # 2. og_presence_absence.tsv
    element_col_keys = sorted(set(element_col_keys))
    n_total_elems = max(1, len(element_col_keys))
    og_cols = [
        "og_id",
        "og_level",
        "cog_category",
        "consensus_description",
        "n_elements",
        "prevalence_pct",
        "classification",
        *element_col_keys,
    ]
    og_rows: list[dict[str, str]] = []
    for og_id, counter in sorted(og_counts.items(), key=lambda kv: (-len(kv[1]), kv[0])):
        n_present = sum(1 for k in element_col_keys if counter.get(k, 0) > 0)
        prev = 100.0 * n_present / n_total_elems
        if prev >= 80.0:
            cls = "CORE_BACKBONE"
        elif prev >= 20.0:
            cls = "SHELL_CARGO"
        else:
            cls = "CLOUD_ACCESSORY"
        meta = og_meta[og_id]
        row_out = {
            "og_id": og_id,
            "og_level": meta["og_level"],
            "cog_category": meta["cog_category"],
            "consensus_description": meta["description"],
            "n_elements": str(n_present),
            "prevalence_pct": f"{prev:.2f}",
            "classification": cls,
        }
        for k in element_col_keys:
            row_out[k] = str(counter.get(k, 0))
        og_rows.append(row_out)

    with open(out_dir / "og_presence_absence.tsv", "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=og_cols, delimiter="\t")
        writer.writeheader()
        writer.writerows(og_rows)

    # 3. MultiQC custom-content YAML/TSV files
    _write_multiqc_files(
        out_dir=out_dir,
        sample_me_counts=sample_me_counts,
        amr_by_sample=amr_by_sample,
        cog_by_sample=cog_by_sample,
        mge_context_counts=mge_context_counts,
    )

    return len(summary_rows), len(og_rows)


def _write_multiqc_files(
    out_dir: Path,
    sample_me_counts: dict[str, Counter[str]],
    amr_by_sample: dict[str, Counter[str]],
    cog_by_sample: dict[str, Counter[str]],
    mge_context_counts: dict[str, Counter[str]],
) -> None:
    # icescreen_summary_mqc.yaml
    yaml_lines = [
        "id: 'advena_icescreen_summary'",
        "section_name: 'advena: ICE/IME Summary'",
        "description: 'Counts of detected Complete ICEs, IMEs, and partial elements per sample.'",
        "plot_type: 'bargraph'",
        "pconfig:",
        "  id: 'advena_icescreen_bar'",
        "  title: 'Detected Mobile Elements per Sample'",
        "  ylab: 'Element count'",
    ]
    if sample_me_counts:
        yaml_lines.append("data:")
        for sample, cnt in sorted(sample_me_counts.items()):
            yaml_lines.append(f"  '{sample}':")
            yaml_lines.append(f"    Complete_ICE: {cnt.get('Complete_ICE', 0)}")
            yaml_lines.append(f"    Complete_IME: {cnt.get('Complete_IME', 0)}")
            yaml_lines.append(f"    Partial_or_Other: {cnt.get('Partial_or_Other', 0)}")
    else:
        yaml_lines.append("data: {}")
    (out_dir / "icescreen_summary_mqc.yaml").write_text(
        "\n".join(yaml_lines) + "\n", encoding="utf-8"
    )

    # amr_card_heatmap_mqc.tsv
    all_amr_genes = sorted({g for cnt in amr_by_sample.values() for g in cnt})
    with open(out_dir / "amr_card_heatmap_mqc.tsv", "w", encoding="utf-8", newline="") as handle:
        handle.write("# id: 'advena_amr_card'\n")
        handle.write("# section_name: 'advena: Dereplicated AMR & Stress Cargo'\n")
        handle.write("# plot_type: 'heatmap'\n")
        if all_amr_genes:
            writer = csv.writer(handle, delimiter="\t")
            writer.writerow(["Sample", *all_amr_genes])
            for sample in sorted(set(sample_me_counts) | set(amr_by_sample)):
                cnt = amr_by_sample.get(sample, Counter())
                writer.writerow([sample, *[cnt.get(g, 0) for g in all_amr_genes]])
        else:
            handle.write("Sample\tNo_AMR_Hits\n")
            for sample in sorted(sample_me_counts):
                handle.write(f"{sample}\t0\n")

    # cog_distribution_mqc.tsv
    all_cogs = sorted({c for cnt in cog_by_sample.values() for c in cnt})
    with open(out_dir / "cog_distribution_mqc.tsv", "w", encoding="utf-8", newline="") as handle:
        handle.write("# id: 'advena_cog_distribution'\n")
        handle.write("# section_name: 'advena: eggNOG COG Functional Categories'\n")
        handle.write("# plot_type: 'bargraph'\n")
        if all_cogs:
            writer = csv.writer(handle, delimiter="\t")
            writer.writerow(["Sample", *all_cogs])
            for sample in sorted(set(sample_me_counts) | set(cog_by_sample)):
                cnt = cog_by_sample.get(sample, Counter())
                writer.writerow([sample, *[cnt.get(c, 0) for c in all_cogs]])
        else:
            handle.write("Sample\tUnassigned\n")
            for sample in sorted(sample_me_counts):
                handle.write(f"{sample}\t0\n")

    # mge_context_mqc.tsv
    ctx_cols = [
        "Plasmid_borne_MEs",
        "Composite_prophage_MEs",
        "IS_flanked_MEs",
        "att_site_resolved_MEs",
    ]
    with open(out_dir / "mge_context_mqc.tsv", "w", encoding="utf-8", newline="") as handle:
        handle.write("# id: 'advena_mge_context'\n")
        handle.write("# section_name: 'advena: Complementary MGE Context'\n")
        handle.write("# plot_type: 'table'\n")
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["Sample", *ctx_cols])
        for sample in sorted(mge_context_counts):
            cnt = mge_context_counts[sample]
            writer.writerow([sample, *[cnt.get(c, 0) for c in ctx_cols]])


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--detected-me", nargs="*", type=Path, default=[])
    parser.add_argument("--contextualized-me", nargs="*", type=Path, default=[])
    parser.add_argument("--cargo-profiles", nargs="*", type=Path, default=[])
    parser.add_argument("--qc-summary", type=Path, default=None)
    parser.add_argument("--clusters", type=Path, default=None)
    parser.add_argument("--outdir", required=True, type=Path)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    n_me, n_og = build_cohort_outputs(
        detected_me_paths=args.detected_me,
        contextualized_me_paths=args.contextualized_me,
        cargo_profile_paths=args.cargo_profiles,
        qc_summary_path=args.qc_summary,
        clusters_path=args.clusters,
        out_dir=args.outdir,
    )
    print(f"Built cohort summary ({n_me} MEs, {n_og} unique OGs)", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
