#!/usr/bin/env python3
"""Merge AMR/CARD, eggNOG OGs, DefenseFinder, and antiSMASH BGCs into master cargo TSV & GFF3."""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path
from urllib.parse import quote

from gbk_utils import format_gff3_attributes, open_text, parse_gff3_attributes

PROFILE_COLUMNS = [
    "sample",
    "me_id",
    "me_category",
    "contig",
    "locus_tag",
    "protein_id",
    "orig_start",
    "orig_end",
    "strand",
    "is_signature_protein",
    "sp_type",
    "product",
    "amr_gene",
    "amr_element_type",
    "amr_drug_class",
    "card_aro",
    "amr_detected_by",
    "eggnog_bacteria_og",
    "eggnog_bacillota_og",
    "eggnog_narrowest_og",
    "cog_category",
    "eggnog_description",
    "cargo_role",
    "defense_system",
    "defense_gene",
    "bgc_type",
]


def _load_tsv_by_key(path: Path | None, key_col: str) -> dict[str, list[dict[str, str]]]:
    grouped: dict[str, list[dict[str, str]]] = {}
    if path is None or not path.exists() or path.stat().st_size == 0:
        return grouped
    with open_text(path, "rt") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            key = (row.get(key_col) or "").strip()
            if key and key != "-":
                grouped.setdefault(key, []).append(row)
    return grouped


def _load_defense_by_locus(
    genes_path: Path | None, systems_path: Path | None
) -> dict[str, tuple[str, str]]:
    """Map locus_tag -> (defense_system_type, defense_gene_name)."""
    sys_type_by_id: dict[str, str] = {}
    if systems_path and systems_path.exists() and systems_path.stat().st_size > 0:
        with open_text(systems_path, "rt") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                sid = (row.get("sys_id") or "").strip()
                stype = (row.get("subtype") or row.get("type") or "").strip()
                if sid and stype:
                    sys_type_by_id[sid] = stype

    by_locus: dict[str, tuple[list[str], list[str]]] = {}
    if genes_path and genes_path.exists() and genes_path.stat().st_size > 0:
        with open_text(genes_path, "rt") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                hit_id = (row.get("hit_id") or "").strip()
                if not hit_id:
                    continue
                locus = hit_id.split()[0]
                sid = (row.get("sys_id") or "").strip()
                stype = sys_type_by_id.get(sid) or (
                    row.get("type") or row.get("subtype") or "Defense"
                )
                gname = (row.get("gene_name") or row.get("hit_gene_ref") or "-").strip()
                systems_list, genes_list = by_locus.setdefault(locus, ([], []))
                if stype and stype not in systems_list:
                    systems_list.append(stype)
                if gname and gname != "-" and gname not in genes_list:
                    genes_list.append(gname)

    return {
        lt: (";".join(sys_l) or "-", ";".join(gene_l) or "-")
        for lt, (sys_l, gene_l) in by_locus.items()
    }


def _load_bgc_by_locus(bgc_path: Path | None) -> tuple[dict[str, list[str]], dict[str, list[str]]]:
    """Return (locus_tag -> [bgc_products], me_id -> [bgc_products])."""
    by_locus: dict[str, list[str]] = {}
    by_me: dict[str, list[str]] = {}
    if bgc_path is None or not bgc_path.exists() or bgc_path.stat().st_size == 0:
        return by_locus, by_me
    with open_text(bgc_path, "rt") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            me_id = (row.get("me_id") or "").strip()
            prods = [p for p in (row.get("bgc_products") or "").split(";") if p and p != "-"]
            if me_id and prods:
                existing = by_me.setdefault(me_id, [])
                for p in prods:
                    if p not in existing:
                        existing.append(p)
            for lt in (row.get("overlapping_locus_tags") or "").split(";"):
                lt = lt.strip()
                if lt and lt != "-":
                    existing_lt = by_locus.setdefault(lt, [])
                    for p in prods:
                        if p not in existing_lt:
                            existing_lt.append(p)
    return by_locus, by_me


def merge_cargo(
    manifest_path: Path,
    icescreen_gff: Path | None,
    amr_path: Path | None,
    eggnog_path: Path | None,
    defense_genes_path: Path | None,
    defense_systems_path: Path | None,
    bgc_path: Path | None,
    out_profile: Path,
    out_gff: Path,
) -> int:
    """Join all cargo layers onto manifest and enrich GFF3."""
    amr_by_locus = _load_tsv_by_key(amr_path, "locus_tag")
    og_by_locus = _load_tsv_by_key(eggnog_path, "locus_tag")
    defense_by_locus = _load_defense_by_locus(defense_genes_path, defense_systems_path)
    bgc_by_locus, bgc_by_me = _load_bgc_by_locus(bgc_path)

    profile_rows: list[dict[str, str]] = []
    annot_by_locus: dict[str, dict[str, str]] = {}

    if manifest_path.exists() and manifest_path.stat().st_size > 0:
        with open_text(manifest_path, "rt") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                lt = (row.get("locus_tag") or "").strip()
                me_id = (row.get("me_id") or "").strip()

                amr_hits = amr_by_locus.get(lt, [])
                amr_gene = ";".join(dict.fromkeys(h["canonical_symbol"] for h in amr_hits)) or "-"
                amr_type = ";".join(dict.fromkeys(h["element_type"] for h in amr_hits)) or "-"
                amr_class = (
                    ";".join(
                        dict.fromkeys(
                            tok
                            for h in amr_hits
                            for tok in h.get("consensus_drug_class", "").split(";")
                            if tok and tok != "-"
                        )
                    )
                    or "-"
                )
                card_aro = (
                    ";".join(
                        dict.fromkeys(
                            h["card_aro_accession"]
                            for h in amr_hits
                            if h.get("card_aro_accession") and h["card_aro_accession"] != "-"
                        )
                    )
                    or "-"
                )
                amr_by = ";".join(dict.fromkeys(h["detected_by"] for h in amr_hits)) or "-"

                og_hits = og_by_locus.get(lt, [])
                og = og_hits[0] if og_hits else {}
                bact_og = og.get("bacteria_og", "-") or "-"
                bac_og = og.get("bacillota_og", "-") or "-"
                narrow_og = og.get("narrowest_og", "-") or "-"
                cog_cat = og.get("cog_category", "-") or "-"
                og_desc = og.get("description", "-") or "-"
                is_sp = row.get("is_signature_protein", "FALSE").upper() == "TRUE"
                cargo_role = og.get("cargo_role") or (
                    "ICE_BACKBONE_SP" if is_sp else "CARGO_UNANNOTATED"
                )

                def_sys, def_gene = defense_by_locus.get(lt, ("-", "-"))
                bgc_prods = bgc_by_locus.get(lt) or []
                bgc_type = ";".join(bgc_prods) if bgc_prods else "-"

                merged = {
                    "sample": row.get("sample", ""),
                    "me_id": me_id,
                    "me_category": row.get("me_category", ""),
                    "contig": row.get("contig", ""),
                    "locus_tag": lt,
                    "protein_id": row.get("protein_id", "-"),
                    "orig_start": row.get("orig_start", "0"),
                    "orig_end": row.get("orig_end", "0"),
                    "strand": row.get("strand", "."),
                    "is_signature_protein": row.get("is_signature_protein", "FALSE"),
                    "sp_type": row.get("sp_type", "-"),
                    "product": row.get("product", "-"),
                    "amr_gene": amr_gene,
                    "amr_element_type": amr_type,
                    "amr_drug_class": amr_class,
                    "card_aro": card_aro,
                    "amr_detected_by": amr_by,
                    "eggnog_bacteria_og": bact_og,
                    "eggnog_bacillota_og": bac_og,
                    "eggnog_narrowest_og": narrow_og,
                    "cog_category": cog_cat,
                    "eggnog_description": og_desc,
                    "cargo_role": cargo_role,
                    "defense_system": def_sys,
                    "defense_gene": def_gene,
                    "bgc_type": bgc_type,
                }
                profile_rows.append(merged)
                if lt:
                    annot_by_locus[lt] = merged

    out_profile.parent.mkdir(parents=True, exist_ok=True)
    with open(out_profile, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=PROFILE_COLUMNS, delimiter="\t")
        writer.writeheader()
        writer.writerows(profile_rows)

    # Emit enriched GFF3 overlaying cargo attributes onto ICEscreen GFF3
    out_gff.parent.mkdir(parents=True, exist_ok=True)
    with open_text(out_gff, "wt") as gff_out:
        if icescreen_gff and icescreen_gff.exists() and icescreen_gff.stat().st_size > 0:
            wrote_header = False
            with open_text(icescreen_gff, "rt") as gff_in:
                for raw_line in gff_in:
                    line = raw_line.rstrip("\r\n")
                    if line.startswith("##gff-version"):
                        wrote_header = True
                    if not line or line.startswith("#"):
                        gff_out.write(line + "\n")
                        continue
                    cols = line.split("\t")
                    if len(cols) < 9:
                        gff_out.write(line + "\n")
                        continue
                    attrs = parse_gff3_attributes(cols[8])
                    lt = attrs.get("locus_tag") or attrs.get("ID") or ""
                    if lt in annot_by_locus:
                        info = annot_by_locus[lt]
                        tag_map = {
                            "amr_gene": info["amr_gene"],
                            "card_aro": info["card_aro"],
                            "eggnog_og": info["eggnog_bacteria_og"],
                            "bacillota_og": info["eggnog_bacillota_og"],
                            "cog_cat": info["cog_category"],
                            "defense_system": info["defense_system"],
                            "bgc_type": info["bgc_type"],
                        }
                        for key, val in tag_map.items():
                            if val and val != "-":
                                attrs[key] = quote(val.replace(";", ","), safe=":,._-|@()")
                        cols[8] = format_gff3_attributes(attrs)
                    else:
                        me_id = attrs.get("ID") or attrs.get("Name") or ""
                        if me_id in bgc_by_me:
                            attrs["bgc_type"] = quote(",".join(bgc_by_me[me_id]), safe=":,._-|@()")
                            cols[8] = format_gff3_attributes(attrs)
                    gff_out.write("\t".join(cols) + "\n")
            if not wrote_header:
                # Input lacked header; valid if empty
                pass
        else:
            gff_out.write("##gff-version 3\n")
            for info in profile_rows:
                attrs = {
                    "ID": info["locus_tag"],
                    "locus_tag": info["locus_tag"],
                    "me_id": info["me_id"],
                }
                for key, col in (
                    ("amr_gene", "amr_gene"),
                    ("card_aro", "card_aro"),
                    ("eggnog_og", "eggnog_bacteria_og"),
                    ("bacillota_og", "eggnog_bacillota_og"),
                    ("cog_cat", "cog_category"),
                    ("defense_system", "defense_system"),
                    ("bgc_type", "bgc_type"),
                ):
                    if info[col] and info[col] != "-":
                        attrs[key] = quote(info[col].replace(";", ","), safe=":,._-|@()")
                gff_out.write(
                    f"{info['contig']}\tadvena\tCDS\t{info['orig_start']}\t{info['orig_end']}\t"
                    f".\t{info['strand']}\t0\t{format_gff3_attributes(attrs)}\n"
                )

    return len(profile_rows)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--icescreen-gff", type=Path, default=None)
    parser.add_argument("--amr", type=Path, default=None)
    parser.add_argument("--eggnog", type=Path, default=None)
    parser.add_argument("--defense-genes", type=Path, default=None)
    parser.add_argument("--defense-systems", type=Path, default=None)
    parser.add_argument("--bgc", type=Path, default=None)
    parser.add_argument("--out-profile", required=True, type=Path)
    parser.add_argument("--out-gff", required=True, type=Path)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    count = merge_cargo(
        manifest_path=args.manifest,
        icescreen_gff=args.icescreen_gff,
        amr_path=args.amr,
        eggnog_path=args.eggnog,
        defense_genes_path=args.defense_genes,
        defense_systems_path=args.defense_systems,
        bgc_path=args.bgc,
        out_profile=args.out_profile,
        out_gff=args.out_gff,
    )
    print(f"Merged cargo profile with {count} CDS row(s)", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
