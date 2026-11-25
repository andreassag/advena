#!/usr/bin/env python3
"""Overlap MOB-suite, geNomad, ISEScan, IntegronFinder, and attL/attR calls with ICE/IMEs."""

from __future__ import annotations

import argparse
import csv
import re
import sys
from pathlib import Path

from gbk_utils import open_text


def _safe_float(val: str | None, default: float = 0.0) -> float:
    if not val:
        return default
    try:
        return float(val.strip())
    except ValueError:
        return default


def _safe_int(val: str | None, default: int = 0) -> int:
    if not val:
        return default
    try:
        return int(float(val.strip()))
    except ValueError:
        return default


def _contig_key(contig: str) -> str:
    return contig.strip().split()[0].split(".")[0]


def _parse_mobsuite(path: Path | None) -> dict[str, tuple[str, str]]:
    """Map contig_key -> (molecule_type, primary_cluster_id)."""
    info: dict[str, tuple[str, str]] = {}
    if path is None or not path.exists() or path.stat().st_size == 0:
        return info
    with open_text(path, "rt") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            cid = (row.get("contig_id") or "").strip()
            mtype = (row.get("molecule_type") or "").strip().lower()
            pcluster = (row.get("primary_cluster_id") or "-").strip()
            if cid:
                info[_contig_key(cid)] = (mtype or "unknown", pcluster)
    return info


def _parse_genomad_plasmids(path: Path | None) -> dict[str, float]:
    """Map contig_key -> plasmid_score."""
    scores: dict[str, float] = {}
    if path is None or not path.exists() or path.stat().st_size == 0:
        return scores
    with open_text(path, "rt") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            seq = (row.get("seq_name") or "").strip()
            score = _safe_float(row.get("plasmid_score"))
            if seq:
                key = _contig_key(seq)
                scores[key] = max(scores.get(key, 0.0), score)
    return scores


def _parse_genomad_viruses(path: Path | None) -> dict[str, list[tuple[int, int, str]]]:
    """Map contig_key -> [(start, end, virus_id)] for provirus/virus regions."""
    regions: dict[str, list[tuple[int, int, str]]] = {}
    if path is None or not path.exists() or path.stat().st_size == 0:
        return regions
    with open_text(path, "rt") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            seq_name = (row.get("seq_name") or "").strip()
            if not seq_name:
                continue
            coords = (row.get("coordinates") or "").strip()
            base_contig = seq_name.split("|provirus")[0]
            start, end = 0, 0
            if coords and "-" in coords:
                parts = coords.split("-", 1)
                start, end = _safe_int(parts[0]), _safe_int(parts[1])
            elif "|provirus_" in seq_name:
                m = re.search(r"provirus_(\d+)_(\d+)", seq_name)
                if m:
                    start, end = int(m.group(1)), int(m.group(2))
            else:
                length = _safe_int(row.get("length"))
                start, end = 1, max(1, length)
            if start > end and end > 0:
                start, end = end, start
            regions.setdefault(_contig_key(base_contig), []).append((start, end, seq_name))
    return regions


def _parse_isescan(path: Path | None) -> dict[str, list[tuple[int, int, str]]]:
    """Map contig_key -> [(start, end, is_family)]."""
    hits: dict[str, list[tuple[int, int, str]]] = {}
    if path is None or not path.exists() or path.stat().st_size == 0:
        return hits
    with open_text(path, "rt") as handle:
        sample_line = handle.readline()
        handle.seek(0)
        delim = "," if sample_line.count(",") > sample_line.count("\t") else "\t"
        for row in csv.DictReader(handle, delimiter=delim):
            seq = (row.get("seqID") or row.get("seqid") or "").strip()
            fam = (row.get("family") or "IS").strip()
            s = _safe_int(row.get("isBegin"))
            e = _safe_int(row.get("isEnd"))
            if seq and e >= s > 0:
                hits.setdefault(_contig_key(seq), []).append((s, e, fam))
    return hits


def _parse_integronfinder(path: Path | None) -> dict[str, list[tuple[int, int, str]]]:
    """Map contig_key -> [(start, end, integron_id)]."""
    integrons: dict[str, dict[str, list[int]]] = {}
    if path is None or not path.exists() or path.stat().st_size == 0:
        return {}
    with open_text(path, "rt") as handle:
        lines = [ln for ln in handle if ln.strip() and not ln.startswith("#")]
    if len(lines) < 2:
        return {}
    reader = csv.DictReader(lines, delimiter="\t")
    for row in reader:
        integ_id = (row.get("ID_integron") or "").strip()
        replicon = (row.get("ID_replicon") or "").strip()
        s = _safe_int(row.get("pos_beg"))
        e = _safe_int(row.get("pos_end"))
        if integ_id and replicon and e > 0:
            key = _contig_key(replicon)
            bounds = integrons.setdefault(key, {}).setdefault(integ_id, [s, e])
            bounds[0] = min(bounds[0], s)
            bounds[1] = max(bounds[1], e)
    return {k: [(b[0], b[1], iid) for iid, b in by_id.items()] for k, by_id in integrons.items()}


def _parse_att_sites(path: Path | None) -> dict[str, dict[str, str]]:
    """Map me_id -> att_site row."""
    by_me: dict[str, dict[str, str]] = {}
    if path is None or not path.exists() or path.stat().st_size == 0:
        return by_me
    with open_text(path, "rt") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            me_id = (row.get("me_id") or "").strip()
            if me_id:
                by_me[me_id] = row
    return by_me


def contextualize_elements(
    detected_me_tsv: Path,
    mobsuite_path: Path | None = None,
    genomad_plasmid_path: Path | None = None,
    genomad_virus_path: Path | None = None,
    isescan_path: Path | None = None,
    integron_path: Path | None = None,
    att_sites_path: Path | None = None,
    is_flank_bp: int = 5000,
    plasmid_score_threshold: float = 0.70,
) -> tuple[list[str], list[dict[str, str]]]:
    """Enrich detected_ME.tsv rows with complementary MGE context columns."""
    extra_cols = [
        "plasmid_borne_warning",
        "mobsuite_molecule_type",
        "mobsuite_cluster_id",
        "genomad_plasmid_score",
        "composite_prophage_ids",
        "flanking_is_families",
        "internal_is_count",
        "integron_ids",
        "att_target_feature",
        "att_repeat_length",
        "att_repeat_seq",
    ]
    if not detected_me_tsv.exists() or detected_me_tsv.stat().st_size == 0:
        return (["ICE_IME_id", "Genome_accession", *extra_cols], [])

    mob = _parse_mobsuite(mobsuite_path)
    g_plas = _parse_genomad_plasmids(genomad_plasmid_path)
    g_vir = _parse_genomad_viruses(genomad_virus_path)
    is_hits = _parse_isescan(isescan_path)
    integrons = _parse_integronfinder(integron_path)
    att_map = _parse_att_sites(att_sites_path)

    with open_text(detected_me_tsv, "rt") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        base_cols = list(reader.fieldnames or [])
        rows = list(reader)

    out_cols = base_cols + [c for c in extra_cols if c not in base_cols]
    enriched: list[dict[str, str]] = []

    for row in rows:
        me_id = (row.get("ICE_IME_id") or "").strip()
        contig = (row.get("Genome_accession") or "").strip()
        ckey = _contig_key(contig)
        me_start = _safe_int(row.get("Start_of_most_upstream_SP"))
        me_end = _safe_int(row.get("Stop_of_most_downstream_SP"))
        if me_start > me_end and me_end > 0:
            me_start, me_end = me_end, me_start

        mob_type, mob_cluster = mob.get(ckey, ("NA", "-"))
        g_score = g_plas.get(ckey)
        is_plasmid = (mob_type == "plasmid") or (
            g_score is not None and g_score >= plasmid_score_threshold
        )

        # Overlapping or nested prophages
        prophages: list[str] = []
        for vs, ve, vid in g_vir.get(ckey, []):
            if me_end > 0 and min(me_end, ve) >= max(me_start, vs):
                prophages.append(vid)

        # Flanking (within is_flank_bp) or internal IS elements
        flanking_is: list[str] = []
        internal_is = 0
        for is_s, is_e, is_fam in is_hits.get(ckey, []):
            if me_end <= 0:
                continue
            if is_s >= me_start and is_e <= me_end:
                internal_is += 1
                if is_fam not in flanking_is:
                    flanking_is.append(is_fam)
            elif (
                min(me_end + is_flank_bp, is_e) >= max(1, me_start - is_flank_bp, is_s)
                and is_fam not in flanking_is
            ):
                flanking_is.append(is_fam)

        # Overlapping integrons
        ov_integrons: list[str] = []
        for int_s, int_e, iid in integrons.get(ckey, []):
            if me_end > 0 and min(me_end, int_e) >= max(me_start, int_s):
                ov_integrons.append(iid)

        att = att_map.get(me_id, {})

        new_row = dict(row)
        new_row.update(
            {
                "plasmid_borne_warning": "TRUE" if is_plasmid else "FALSE",
                "mobsuite_molecule_type": mob_type,
                "mobsuite_cluster_id": mob_cluster,
                "genomad_plasmid_score": f"{g_score:.4f}" if g_score is not None else "NA",
                "composite_prophage_ids": ";".join(dict.fromkeys(prophages)) or "-",
                "flanking_is_families": ";".join(flanking_is) or "-",
                "internal_is_count": str(internal_is),
                "integron_ids": ";".join(dict.fromkeys(ov_integrons)) or "-",
                "att_target_feature": att.get("att_target_feature", "-") or "-",
                "att_repeat_length": att.get("att_repeat_length", "0") or "0",
                "att_repeat_seq": att.get("att_repeat_seq", "-") or "-",
            }
        )
        enriched.append(new_row)

    return out_cols, enriched


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--detected-me", required=True, type=Path)
    parser.add_argument("--mobsuite", type=Path, default=None)
    parser.add_argument("--genomad-plasmid", type=Path, default=None)
    parser.add_argument("--genomad-virus", type=Path, default=None)
    parser.add_argument("--isescan", type=Path, default=None)
    parser.add_argument("--integronfinder", type=Path, default=None)
    parser.add_argument("--att-sites", type=Path, default=None)
    parser.add_argument("--is-flank-bp", type=int, default=5000)
    parser.add_argument("--plasmid-threshold", type=float, default=0.70)
    parser.add_argument("-o", "--output", required=True, type=Path)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    cols, rows = contextualize_elements(
        detected_me_tsv=args.detected_me,
        mobsuite_path=args.mobsuite,
        genomad_plasmid_path=args.genomad_plasmid,
        genomad_virus_path=args.genomad_virus,
        isescan_path=args.isescan,
        integron_path=args.integronfinder,
        att_sites_path=args.att_sites,
        is_flank_bp=args.is_flank_bp,
        plasmid_score_threshold=args.plasmid_threshold,
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with open(args.output, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=cols, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    print(f"Contextualized {len(rows)} mobile element(s)", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
