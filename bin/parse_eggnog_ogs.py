#!/usr/bin/env python3
"""Parse eggNOG-mapper (.emapper.annotations) into hierarchical OG tables for mobile cargo."""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

from gbk_utils import open_text

OUTPUT_COLUMNS = [
    "sample",
    "me_id",
    "contig",
    "locus_tag",
    "is_signature_protein",
    "cargo_role",
    "root_cog",
    "bacteria_og",
    "bacillota_og",
    "narrowest_og",
    "cog_category",
    "preferred_name",
    "description",
    "kegg_ko",
    "kegg_pathway",
    "pfams",
    "gos",
    "seed_ortholog",
    "evalue",
    "score",
]


def parse_og_hierarchy(og_string: str) -> dict[str, str]:
    """Extract @1 (root), @2 (Bacteria), @1239 (Bacillota/Firmicutes), and narrowest OG."""
    result = {
        "root_cog": "-",
        "bacteria_og": "-",
        "bacillota_og": "-",
        "narrowest_og": "-",
    }
    if not og_string or og_string.strip() in {"-", "NA"}:
        return result

    entries = [tok.strip() for tok in og_string.split(",") if tok.strip()]
    ordered_ogs: list[str] = []
    for entry in entries:
        if "@" not in entry:
            continue
        og_id, tax_part = entry.split("@", 1)
        tax_id = tax_part.split("|", 1)[0].strip()
        og_token = f"{og_id}@{tax_id}"
        ordered_ogs.append(og_token)
        if tax_id == "1":
            result["root_cog"] = og_id
        elif tax_id == "2":
            result["bacteria_og"] = og_token
        elif tax_id == "1239":
            result["bacillota_og"] = og_token

    if ordered_ogs:
        result["narrowest_og"] = ordered_ogs[-1]
    return result


def classify_cargo_role(is_sp: bool, cog_category: str, has_og: bool) -> str:
    """Classify CDS into ICE_BACKBONE_SP, CARGO_CHARACTERIZED_OG, or CARGO_HYPOTHETICAL_OG."""
    if is_sp:
        return "ICE_BACKBONE_SP"
    cleaned_cog = cog_category.strip().upper()
    if not has_og or cleaned_cog in {"", "-", "S"}:
        return "CARGO_HYPOTHETICAL_OG"
    return "CARGO_CHARACTERIZED_OG"


def parse_emapper_annotations(
    emapper_path: Path | None,
    manifest_path: Path | None,
    sample: str,
) -> list[dict[str, str]]:
    """Join .emapper.annotations onto *_ME_cds_manifest.tsv."""
    manifest: dict[str, dict[str, str]] = {}
    manifest_order: list[str] = []
    if manifest_path and manifest_path.exists() and manifest_path.stat().st_size > 0:
        with open_text(manifest_path, "rt") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                lt = (row.get("locus_tag") or "").strip()
                if lt:
                    manifest[lt] = row
                    manifest_order.append(lt)

    emapper_by_query: dict[str, dict[str, str]] = {}
    if emapper_path and emapper_path.exists() and emapper_path.stat().st_size > 0:
        with open_text(emapper_path, "rt") as handle:
            header: list[str] | None = None
            for raw_line in handle:
                line = raw_line.rstrip("\r\n")
                if not line:
                    continue
                if line.startswith("##"):
                    continue
                if line.startswith("#query"):
                    header = [col.lstrip("#").strip() for col in line.split("\t")]
                    continue
                if line.startswith("#"):
                    continue
                if header is None:
                    # Fallback if file has header without leading '#'
                    header = [col.strip() for col in line.split("\t")]
                    if header and header[0] == "query":
                        continue
                fields = line.split("\t")
                row = {
                    header[i]: (fields[i].strip() if i < len(fields) else "-")
                    for i in range(len(header))
                }
                q = row.get("query", "").split()[0]
                if q:
                    emapper_by_query[q] = row

    all_queries = list(manifest_order)
    for q in emapper_by_query:
        if q not in manifest:
            all_queries.append(q)

    rows: list[dict[str, str]] = []
    for locus_tag in all_queries:
        m_info = manifest.get(locus_tag, {})
        e_info = emapper_by_query.get(locus_tag, {})
        og_parsed = parse_og_hierarchy(e_info.get("eggNOG_OGs", ""))
        cog_cat = e_info.get("COG_category", "-") or "-"
        is_sp = m_info.get("is_signature_protein", "FALSE").upper() == "TRUE"
        has_og = any(v != "-" for v in og_parsed.values())
        role = classify_cargo_role(is_sp=is_sp, cog_category=cog_cat, has_og=has_og)

        rows.append(
            {
                "sample": sample,
                "me_id": m_info.get("me_id", "UNASSIGNED"),
                "contig": m_info.get("contig", "UNKNOWN"),
                "locus_tag": locus_tag,
                "is_signature_protein": "TRUE" if is_sp else "FALSE",
                "cargo_role": role,
                "root_cog": og_parsed["root_cog"],
                "bacteria_og": og_parsed["bacteria_og"],
                "bacillota_og": og_parsed["bacillota_og"],
                "narrowest_og": og_parsed["narrowest_og"],
                "cog_category": cog_cat,
                "preferred_name": e_info.get("Preferred_name", "-") or "-",
                "description": e_info.get("Description", "-") or "-",
                "kegg_ko": e_info.get("KEGG_ko", "-") or "-",
                "kegg_pathway": e_info.get("KEGG_Pathway", "-") or "-",
                "pfams": e_info.get("PFAMs", "-") or "-",
                "gos": e_info.get("GOs", "-") or "-",
                "seed_ortholog": e_info.get("seed_ortholog", "-") or "-",
                "evalue": e_info.get("evalue", "-") or "-",
                "score": e_info.get("score", "-") or "-",
            }
        )
    return rows


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sample", required=True, help="Sample identifier")
    parser.add_argument("--emapper", type=Path, default=None, help=".emapper.annotations file")
    parser.add_argument("--manifest", type=Path, default=None, help="*_ME_cds_manifest.tsv")
    parser.add_argument("-o", "--output", required=True, type=Path, help="Output cargo OGs TSV")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    rows = parse_emapper_annotations(args.emapper, args.manifest, args.sample)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with open(args.output, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=OUTPUT_COLUMNS, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    print(f"Parsed {len(rows)} eggNOG OG row(s) for {args.sample}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
