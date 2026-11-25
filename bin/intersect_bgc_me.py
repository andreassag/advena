#!/usr/bin/env python3
"""Intersect antiSMASH-lite secondary metabolite BGC regions with detected ICE/IME intervals."""

from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

from gbk_utils import open_text, parse_genbank

OUTPUT_COLUMNS = [
    "sample",
    "me_id",
    "contig",
    "me_start",
    "me_end",
    "bgc_region_id",
    "bgc_products",
    "bgc_start",
    "bgc_end",
    "overlap_bp",
    "me_overlap_frac",
    "bgc_overlap_frac",
    "overlapping_locus_tags",
]


def _extract_bgcs_from_json(json_path: Path) -> list[dict]:
    """Extract BGC region intervals from an antiSMASH JSON output file."""
    bgcs: list[dict] = []
    if not json_path.exists() or json_path.stat().st_size == 0:
        return bgcs
    data = json.loads(json_path.read_text(encoding="utf-8"))
    for rec_idx, record in enumerate(data.get("records", []), start=1):
        contig = str(record.get("id") or record.get("name") or f"record_{rec_idx}")
        areas = record.get("areas", [])
        for area_idx, area in enumerate(areas, start=1):
            start = int(area.get("start", 0)) + 1  # antiSMASH JSON areas are 0-based
            end = int(area.get("end", 0))
            products = area.get("products", [])
            if isinstance(products, str):
                products = [products]
            bgcs.append(
                {
                    "contig": contig,
                    "region_id": f"{contig}.region{area_idx:03d}",
                    "products": sorted({str(p) for p in products if p}),
                    "start": start,
                    "end": end,
                }
            )
        # Fallback: scan features for 'region' entries if areas absent
        if not areas:
            region_idx = 0
            for feat in record.get("features", []):
                if feat.get("type") != "region":
                    continue
                region_idx += 1
                loc = str(feat.get("location", ""))
                nums = [int(x) for x in __import__("re").findall(r"\d+", loc)]
                if len(nums) >= 2:
                    quals = feat.get("qualifiers", {})
                    prods = quals.get("product", ["unknown"])
                    bgcs.append(
                        {
                            "contig": contig,
                            "region_id": f"{contig}.region{region_idx:03d}",
                            "products": sorted({str(p) for p in prods}),
                            "start": min(nums) + 1,
                            "end": max(nums),
                        }
                    )
    return bgcs


def _extract_bgcs_from_gbk(gbk_paths: list[Path]) -> list[dict]:
    """Extract BGC region intervals from antiSMASH annotated GenBank files."""
    bgcs: list[dict] = []
    for gbk_path in gbk_paths:
        if not gbk_path.exists() or gbk_path.stat().st_size == 0:
            continue
        for rec in parse_genbank(gbk_path):
            region_idx = 0
            for feat in rec.features:
                if feat.feature_type not in {"region", "cand_cluster", "protocluster"}:
                    continue
                if feat.feature_type != "region":
                    continue
                region_idx += 1
                prods = feat.qualifiers.get("product", ["unknown"])
                bgcs.append(
                    {
                        "contig": rec.record_id,
                        "region_id": f"{rec.record_id}.region{region_idx:03d}",
                        "products": sorted({p for p in prods if p}),
                        "start": feat.start,
                        "end": feat.end,
                    }
                )
    return bgcs


def intersect_bgcs(
    sample: str,
    detected_me_tsv: Path,
    manifest_tsv: Path | None,
    antismash_json: Path | None,
    antismash_gbks: list[Path],
) -> list[dict[str, str]]:
    """Compute coordinate overlaps between antiSMASH BGCs and detected ICE/IME intervals."""
    bgcs = _extract_bgcs_from_json(antismash_json) if antismash_json else []
    if not bgcs and antismash_gbks:
        bgcs = _extract_bgcs_from_gbk(antismash_gbks)
    if not bgcs or not detected_me_tsv.exists() or detected_me_tsv.stat().st_size == 0:
        return []

    cds_by_contig: dict[str, list[tuple[int, int, str]]] = {}
    if manifest_tsv and manifest_tsv.exists() and manifest_tsv.stat().st_size > 0:
        with open_text(manifest_tsv, "rt") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                contig = (row.get("contig") or "").strip()
                lt = (row.get("locus_tag") or "").strip()
                try:
                    s, e = int(row["orig_start"]), int(row["orig_end"])
                except (KeyError, ValueError):
                    continue
                cds_by_contig.setdefault(contig, []).append((s, e, lt))

    overlaps: list[dict[str, str]] = []
    with open_text(detected_me_tsv, "rt") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            me_id = (row.get("ICE_IME_id") or "").strip()
            contig = (row.get("Genome_accession") or "").strip()
            if not me_id or me_id == "-":
                continue
            try:
                me_start = int(float(row["Start_of_most_upstream_SP"]))
                me_end = int(float(row["Stop_of_most_downstream_SP"]))
            except (KeyError, ValueError):
                continue
            if me_start > me_end:
                me_start, me_end = me_end, me_start

            me_len = max(1, me_end - me_start + 1)
            for bgc in bgcs:
                if bgc["contig"].split(".")[0] != contig.split(".")[0]:
                    continue
                ov = max(0, min(me_end, bgc["end"]) - max(me_start, bgc["start"]) + 1)
                if ov <= 0:
                    continue
                bgc_len = max(1, bgc["end"] - bgc["start"] + 1)
                ov_start = max(me_start, bgc["start"])
                ov_end = min(me_end, bgc["end"])
                loci = [
                    lt
                    for s, e, lt in cds_by_contig.get(contig, [])
                    if min(e, ov_end) >= max(s, ov_start)
                ]
                overlaps.append(
                    {
                        "sample": sample,
                        "me_id": me_id,
                        "contig": contig,
                        "me_start": str(me_start),
                        "me_end": str(me_end),
                        "bgc_region_id": bgc["region_id"],
                        "bgc_products": ";".join(bgc["products"]) or "unknown",
                        "bgc_start": str(bgc["start"]),
                        "bgc_end": str(bgc["end"]),
                        "overlap_bp": str(ov),
                        "me_overlap_frac": f"{ov / me_len:.4f}",
                        "bgc_overlap_frac": f"{ov / bgc_len:.4f}",
                        "overlapping_locus_tags": ";".join(loci) if loci else "-",
                    }
                )
    return overlaps


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--detected-me", required=True, type=Path)
    parser.add_argument("--manifest", type=Path, default=None)
    parser.add_argument("--antismash-json", type=Path, default=None)
    parser.add_argument("--antismash-gbk", nargs="*", type=Path, default=[])
    parser.add_argument("-o", "--output", required=True, type=Path)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    rows = intersect_bgcs(
        sample=args.sample,
        detected_me_tsv=args.detected_me,
        manifest_tsv=args.manifest,
        antismash_json=args.antismash_json,
        antismash_gbks=args.antismash_gbk,
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with open(args.output, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=OUTPUT_COLUMNS, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    print(f"Found {len(rows)} BGC-ME overlap(s) for {args.sample}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
