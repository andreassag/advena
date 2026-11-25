#!/usr/bin/env python3
"""Slice detected ICE/IME intervals from annotated GenBank into BED, FNA, FAA, GFF3, and GBK."""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

from gbk_utils import (
    GBFeature,
    GBRecord,
    extract_cds_protein,
    format_genbank_record,
    open_text,
    parse_genbank,
    wrap_fasta,
)


def _match_record(records_by_id: dict[str, GBRecord], accession: str) -> GBRecord | None:
    if accession in records_by_id:
        return records_by_id[accession]
    base = accession.split(".")[0]
    for rec in records_by_id.values():
        if rec.accession == accession or rec.accession == base or rec.locus == base:
            return rec
    return None


def _parse_sp_map(sp_tsv: Path | None) -> dict[str, str]:
    """Map CDS_locus_tag -> CDS_Protein_type for detected signature proteins."""
    sp_map: dict[str, str] = {}
    if sp_tsv is None or not sp_tsv.exists() or sp_tsv.stat().st_size == 0:
        return sp_map
    with open_text(sp_tsv, "rt") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            locus = (row.get("CDS_locus_tag") or "").strip()
            ptype = (row.get("CDS_Protein_type") or "").strip()
            if locus and locus != "-":
                sp_map[locus] = ptype or "SP"
    return sp_map


def _parse_detected_me(me_tsv: Path) -> list[dict[str, str]]:
    """Read *_detected_ME.tsv rows, returning empty list if header-only or empty."""
    if not me_tsv.exists() or me_tsv.stat().st_size == 0:
        return []
    rows: list[dict[str, str]] = []
    with open_text(me_tsv, "rt") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            me_id = (row.get("ICE_IME_id") or "").strip()
            if not me_id or me_id == "-":
                continue
            rows.append(row)
    return rows


def extract_me_sequences(
    sample: str,
    prefix: str,
    detected_me_tsv: Path,
    sp_with_me_tsv: Path | None,
    genbank_path: Path,
    out_dir: Path,
    flank_bp: int = 0,
    scope: str = "elements",
) -> int:
    """Extract ICE/IME sequences, CDS proteomes, GFF3, manifest, and per-ME GenBank slices."""
    if flank_bp < 0:
        raise ValueError("flank_bp must be >= 0")

    records = parse_genbank(genbank_path)
    records_by_id = {rec.record_id: rec for rec in records}
    sp_map = _parse_sp_map(sp_with_me_tsv)
    me_rows = _parse_detected_me(detected_me_tsv)

    out_dir.mkdir(parents=True, exist_ok=True)
    slices_dir = out_dir / "slices"
    slices_dir.mkdir(parents=True, exist_ok=True)

    bed_path = out_dir / f"{prefix}_ME_regions.bed"
    fna_path = out_dir / f"{prefix}_ME_regions.fna.gz"
    cargo_fna_path = out_dir / f"{prefix}_cargo_query.fna.gz"
    faa_path = out_dir / f"{prefix}_ME_cds.faa.gz"
    gff_path = out_dir / f"{prefix}_ME_cds.gff3.gz"
    manifest_path = out_dir / f"{prefix}_ME_cds_manifest.tsv"
    count_path = out_dir / f"{prefix}_me_count.txt"

    bed_lines: list[str] = []
    fna_entries: list[str] = []
    cargo_fna_entries: list[str] = []
    faa_entries: list[str] = []
    gff_lines: list[str] = ["##gff-version 3"]
    manifest_rows: list[dict[str, str]] = []

    seen_locus_tags: set[str] = set()
    extracted_elements = 0

    def _assign_locus_tag(rec: GBRecord, feat: GBFeature, idx: int) -> str:
        raw = (
            feat.get_first("locus_tag")
            or feat.get_first("protein_id")
            or f"{rec.record_id}_CDS_{idx}"
        )
        tag = raw.strip().replace(" ", "_")
        if tag in seen_locus_tags:
            suffix = 2
            while f"{tag}_{suffix}" in seen_locus_tags:
                suffix += 1
            tag = f"{tag}_{suffix}"
        seen_locus_tags.add(tag)
        return tag

    # Map CDS -> owning ME for genome-scope tagging
    cds_to_me: dict[tuple[str, int, int], tuple[str, str]] = {}

    for me_row in me_rows:
        me_id = me_row["ICE_IME_id"].strip()
        contig = (me_row.get("Genome_accession") or "").strip()
        category = (me_row.get("Category_of_element") or "Unknown").strip()
        try:
            raw_start = int(float(me_row["Start_of_most_upstream_SP"]))
            raw_end = int(float(me_row["Stop_of_most_downstream_SP"]))
        except (KeyError, ValueError):
            continue
        if raw_start > raw_end:
            raw_start, raw_end = raw_end, raw_start

        rec = _match_record(records_by_id, contig)
        if rec is None or not rec.sequence:
            print(
                f"WARNING: Contig '{contig}' for {me_id} not found in {genbank_path}; skipping.",
                file=sys.stderr,
            )
            continue

        contig_len = len(rec.sequence)
        slice_start = max(1, raw_start - flank_bp)
        slice_end = min(contig_len, raw_end + flank_bp)
        slice_seq = rec.sequence[slice_start - 1 : slice_end]
        slice_name = f"{sample}__{me_id}"
        extracted_elements += 1

        # BED6 (0-based half-open)
        bed_lines.append(
            f"{rec.record_id}\t{slice_start - 1}\t{slice_end}\t{me_id}\t0\t+\t{category}"
        )

        # Element FNA
        fna_entries.append(
            f">{slice_name} contig={rec.record_id} start={slice_start} end={slice_end} "
            f"category={category}\n{wrap_fasta(slice_seq)}\n"
        )

        # Collect features within slice for per-element GBK and CDS outputs
        shifted_features: list[GBFeature] = []
        cds_idx = 0
        for feat in rec.features:
            if feat.feature_type == "source":
                continue
            if feat.end < slice_start or feat.start > slice_end:
                continue

            shifted_start = max(1, feat.start - slice_start + 1)
            shifted_end = min(len(slice_seq), feat.end - slice_start + 1)
            shifted_quals = {k: list(v) for k, v in feat.qualifiers.items()}

            if feat.feature_type == "CDS":
                cds_idx += 1
                orig_locus = feat.get_first("locus_tag") or feat.get_first("protein_id")
                is_sp = orig_locus in sp_map
                sp_type = sp_map.get(orig_locus, "-")
                if is_sp:
                    shifted_quals.setdefault("note", []).append(f"ICEscreen_SP:{sp_type}")
                cds_to_me[(rec.record_id, feat.start, feat.end)] = (me_id, category)

                if scope == "elements":
                    locus_tag = _assign_locus_tag(rec, feat, cds_idx)
                    protein_id = feat.get_first("protein_id") or "-"
                    product = feat.get_first("product") or "hypothetical protein"
                    protein = extract_cds_protein(rec, feat)
                    if protein:
                        faa_entries.append(
                            f">{locus_tag} me_id={me_id} contig={rec.record_id} "
                            f"start={feat.start} end={feat.end} strand={feat.strand} "
                            f"{product}\n{wrap_fasta(protein)}\n"
                        )
                    gff_lines.append(
                        f"{slice_name}\tadvena\tCDS\t{shifted_start}\t{shifted_end}\t.\t"
                        f"{feat.strand}\t0\tID={locus_tag};Name={locus_tag};"
                        f"locus_tag={locus_tag};me_id={me_id}"
                    )
                    manifest_rows.append(
                        {
                            "sample": sample,
                            "me_id": me_id,
                            "me_category": category,
                            "contig": rec.record_id,
                            "locus_tag": locus_tag,
                            "protein_id": protein_id,
                            "orig_start": str(feat.start),
                            "orig_end": str(feat.end),
                            "strand": feat.strand,
                            "is_signature_protein": "TRUE" if is_sp else "FALSE",
                            "sp_type": sp_type,
                            "product": product,
                        }
                    )

            shifted_features.append(
                GBFeature(
                    feature_type=feat.feature_type,
                    location=f"{shifted_start}..{shifted_end}",
                    start=shifted_start,
                    end=shifted_end,
                    strand=feat.strand,
                    parts=[(shifted_start, shifted_end)],
                    qualifiers=shifted_quals,
                )
            )

        gbk_text = format_genbank_record(
            locus_name=slice_name,
            definition=(
                f"{sample} {me_id} ({category}) {rec.record_id}:{slice_start}..{slice_end}"
            ),
            organism=rec.organism or "Bacteria",
            sequence=slice_seq,
            features=shifted_features,
            topology="linear",
        )
        with open_text(slices_dir / f"{slice_name}.gbk.gz", "wt") as gbk_out:
            gbk_out.write(gbk_text)

    if scope == "genome":
        # Whole-genome CDS scope: include all CDS, tagging those inside an ME,
        # and emit full contig nucleotide FASTA matching the contig-coordinate GFF3
        seen_locus_tags.clear()
        for rec in records:
            if rec.sequence:
                cargo_fna_entries.append(f">{rec.record_id}\n{wrap_fasta(rec.sequence)}\n")
            cds_idx = 0
            for feat in rec.features:
                if feat.feature_type != "CDS":
                    continue
                cds_idx += 1
                orig_locus = feat.get_first("locus_tag") or feat.get_first("protein_id")
                locus_tag = _assign_locus_tag(rec, feat, cds_idx)
                protein_id = feat.get_first("protein_id") or "-"
                product = feat.get_first("product") or "hypothetical protein"
                is_sp = orig_locus in sp_map
                sp_type = sp_map.get(orig_locus, "-")
                me_id, category = cds_to_me.get(
                    (rec.record_id, feat.start, feat.end), ("GENOME_BG", "Host_background")
                )
                protein = extract_cds_protein(rec, feat)
                if protein:
                    faa_entries.append(
                        f">{locus_tag} me_id={me_id} contig={rec.record_id} "
                        f"start={feat.start} end={feat.end} strand={feat.strand} "
                        f"{product}\n{wrap_fasta(protein)}\n"
                    )
                gff_lines.append(
                    f"{rec.record_id}\tadvena\tCDS\t{feat.start}\t{feat.end}\t.\t"
                    f"{feat.strand}\t0\tID={locus_tag};Name={locus_tag};"
                    f"locus_tag={locus_tag};me_id={me_id}"
                )
                manifest_rows.append(
                    {
                        "sample": sample,
                        "me_id": me_id,
                        "me_category": category,
                        "contig": rec.record_id,
                        "locus_tag": locus_tag,
                        "protein_id": protein_id,
                        "orig_start": str(feat.start),
                        "orig_end": str(feat.end),
                        "strand": feat.strand,
                        "is_signature_protein": "TRUE" if is_sp else "FALSE",
                        "sp_type": sp_type,
                        "product": product,
                    }
                )
    else:
        cargo_fna_entries = fna_entries

    bed_path.write_text("\n".join(bed_lines) + ("\n" if bed_lines else ""), encoding="utf-8")
    with open_text(fna_path, "wt") as handle:
        handle.write("".join(fna_entries))
    with open_text(cargo_fna_path, "wt") as handle:
        handle.write("".join(cargo_fna_entries))
    with open_text(faa_path, "wt") as handle:
        handle.write("".join(faa_entries))
    with open_text(gff_path, "wt") as handle:
        handle.write("\n".join(gff_lines) + "\n")

    manifest_cols = [
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
    ]
    with open(manifest_path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=manifest_cols, delimiter="\t")
        writer.writeheader()
        writer.writerows(manifest_rows)

    count_path.write_text(f"{extracted_elements}\n", encoding="utf-8")
    return extracted_elements


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sample", required=True, help="Sample identifier")
    parser.add_argument("--prefix", required=True, help="Output file prefix")
    parser.add_argument("--detected-me", required=True, type=Path, help="*_detected_ME.tsv")
    parser.add_argument("--sp-with-me", type=Path, default=None, help="*_detected_SP_withMEIds.tsv")
    parser.add_argument("--genbank", required=True, type=Path, help="*_icescreen.gb.gz or .gbff")
    parser.add_argument("--outdir", required=True, type=Path, help="Output directory")
    parser.add_argument("--flank-bp", type=int, default=0, help="Flanking bp around ME bounds")
    parser.add_argument(
        "--scope",
        choices=("elements", "genome"),
        default="elements",
        help="Emit CDS for extracted elements only or whole genome",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    count = extract_me_sequences(
        sample=args.sample,
        prefix=args.prefix,
        detected_me_tsv=args.detected_me,
        sp_with_me_tsv=args.sp_with_me,
        genbank_path=args.genbank,
        out_dir=args.outdir,
        flank_bp=args.flank_bp,
        scope=args.scope,
    )
    print(f"Extracted {count} mobile element(s) for {args.sample}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
