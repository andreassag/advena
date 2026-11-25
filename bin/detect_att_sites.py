#!/usr/bin/env python3
"""Detect attL/attR direct repeats and tRNA/housekeeping targets flanking ICE/IME integrases."""

from __future__ import annotations

import argparse
import csv
import shutil
import subprocess
import sys
import tempfile
from dataclasses import dataclass
from pathlib import Path

from gbk_utils import open_text, parse_genbank, parse_gff3_attributes, wrap_fasta

HOUSEKEEPING_TARGETS = {"rpsi", "rlmh", "guaa", "ruma", "fda", "guab", "lyss", "ssra"}

OUTPUT_COLUMNS = [
    "sample",
    "me_id",
    "contig",
    "integrase_side",
    "attl_start",
    "attl_end",
    "attr_start",
    "attr_end",
    "att_repeat_length",
    "att_identity_pct",
    "att_repeat_seq",
    "att_target_feature",
    "att_target_distance_bp",
    "detection_method",
]


@dataclass
class RepeatCandidate:
    attl_start: int
    attl_end: int
    attr_start: int
    attr_end: int
    length: int
    identity_pct: float
    seq: str
    method: str


@dataclass
class TargetFeature:
    contig: str
    feature_type: str
    name: str
    start: int
    end: int
    strand: str


def _parse_targets_from_gff(gff_path: Path | None) -> dict[str, list[TargetFeature]]:
    targets: dict[str, list[TargetFeature]] = {}
    if gff_path is None or not gff_path.exists() or gff_path.stat().st_size == 0:
        return targets
    with open_text(gff_path, "rt") as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            cols = line.rstrip("\r\n").split("\t")
            if len(cols) < 9:
                continue
            contig, _src, ftype, start_s, end_s, _score, strand, _phase, attr_s = cols
            attrs = parse_gff3_attributes(attr_s)
            gene = attrs.get("gene") or attrs.get("Name") or attrs.get("locus_tag") or ""
            product = attrs.get("product") or ""
            ftype_low = ftype.lower()
            is_rna = ftype_low in {"trna", "tmrna", "rrna"}
            is_hk = gene.lower() in HOUSEKEEPING_TARGETS or any(
                hk in product.lower() for hk in ("tmrna", "trna-")
            )
            if not (is_rna or is_hk):
                continue
            try:
                s, e = int(start_s), int(end_s)
            except ValueError:
                continue
            label = f"{ftype}:{gene or product or 'target'}"
            targets.setdefault(contig, []).append(
                TargetFeature(
                    contig=contig,
                    feature_type=ftype,
                    name=label,
                    start=min(s, e),
                    end=max(s, e),
                    strand=strand,
                )
            )
    return targets


def _parse_targets_from_gbk(records_by_id: dict) -> dict[str, list[TargetFeature]]:
    targets: dict[str, list[TargetFeature]] = {}
    for contig, rec in records_by_id.items():
        for feat in rec.features:
            ftype = feat.feature_type
            gene = feat.get_first("gene") or feat.get_first("locus_tag") or ""
            product = feat.get_first("product") or ""
            is_rna = ftype.lower() in {"trna", "tmrna", "rrna"}
            is_hk = gene.lower() in HOUSEKEEPING_TARGETS
            if not (is_rna or is_hk):
                continue
            label = f"{ftype}:{gene or product or 'target'}"
            targets.setdefault(contig, []).append(
                TargetFeature(
                    contig=contig,
                    feature_type=ftype,
                    name=label,
                    start=feat.start,
                    end=feat.end,
                    strand=feat.strand,
                )
            )
    return targets


def _is_low_complexity(seq: str) -> bool:
    """Reject homopolymer / dinucleotide low-complexity tracts as att repeat artifacts."""
    if not seq:
        return True
    if len(set(seq)) < 3:
        return True
    max_freq = max(seq.count(base) for base in "ACGT") / len(seq)
    return max_freq >= 0.75


def _merge_targets(
    primary: dict[str, list[TargetFeature]],
    secondary: dict[str, list[TargetFeature]],
) -> dict[str, list[TargetFeature]]:
    """Combine targets from GFF and GenBank, deduplicating by (contig, start, end)."""
    merged: dict[str, list[TargetFeature]] = {k: list(v) for k, v in primary.items()}
    seen: set[tuple[str, int, int]] = {
        (contig, f.start, f.end) for contig, feats in merged.items() for f in feats
    }
    for contig, feats in secondary.items():
        for feat in feats:
            key = (contig, feat.start, feat.end)
            if key not in seen:
                seen.add(key)
                merged.setdefault(contig, []).append(feat)
    return merged


def find_direct_repeats_kmer(
    left_seq: str,
    left_offset: int,
    right_seq: str,
    right_offset: int,
    min_len: int = 8,
    max_len: int = 65,
    min_identity: float = 85.0,
    seed_k: int = 7,
) -> list[RepeatCandidate]:
    """Pure-Python seed-and-extend direct repeat finder (same orientation, 8-65 bp)."""
    l_up = left_seq.upper()
    r_up = right_seq.upper()
    if len(l_up) < min_len or len(r_up) < min_len:
        return []

    r_index: dict[str, list[int]] = {}
    for pos in range(0, len(r_up) - seed_k + 1):
        kmer = r_up[pos : pos + seed_k]
        if "N" in kmer:
            continue
        r_index.setdefault(kmer, []).append(pos)

    seen_spans: set[tuple[int, int, int, int]] = set()
    candidates: list[RepeatCandidate] = []

    for l_pos in range(0, len(l_up) - seed_k + 1):
        kmer = l_up[l_pos : l_pos + seed_k]
        r_positions = r_index.get(kmer)
        if not r_positions:
            continue
        for r_pos in r_positions:
            # Skip redundant interior seeds that extend identically from the left
            if l_pos > 0 and r_pos > 0 and l_up[l_pos - 1] == r_up[r_pos - 1]:
                continue

            # Extend right allowing mismatches while keeping identity >= min_identity
            best_len = seed_k
            best_matches = seed_k
            matches = seed_k
            cur_len = seed_k
            max_possible = min(max_len, len(l_up) - l_pos, len(r_up) - r_pos)
            mismatches_in_row = 0

            while cur_len < max_possible:
                lc = l_up[l_pos + cur_len]
                rc = r_up[r_pos + cur_len]
                cur_len += 1
                if lc == rc and lc != "N":
                    matches += 1
                    mismatches_in_row = 0
                else:
                    mismatches_in_row += 1
                    if mismatches_in_row > 2:
                        break
                ident = 100.0 * matches / cur_len
                if cur_len >= min_len and ident >= min_identity and lc == rc:
                    best_len = cur_len
                    best_matches = matches

            if best_len < min_len:
                continue
            ident = 100.0 * best_matches / best_len
            if ident < min_identity:
                continue

            rep_seq = l_up[l_pos : l_pos + best_len]
            if _is_low_complexity(rep_seq):
                continue

            attl_s = left_offset + l_pos
            attl_e = attl_s + best_len - 1
            attr_s = right_offset + r_pos
            attr_e = attr_s + best_len - 1
            key = (attl_s, attl_e, attr_s, attr_e)
            if key in seen_spans:
                continue
            seen_spans.add(key)
            candidates.append(
                RepeatCandidate(
                    attl_start=attl_s,
                    attl_end=attl_e,
                    attr_start=attr_s,
                    attr_end=attr_e,
                    length=best_len,
                    identity_pct=round(ident, 2),
                    seq=rep_seq,
                    method="kmer_direct_repeat",
                )
            )
    return candidates


def find_direct_repeats_blastn(
    left_seq: str,
    left_offset: int,
    right_seq: str,
    right_offset: int,
    min_len: int = 8,
    max_len: int = 65,
    min_identity: float = 85.0,
) -> list[RepeatCandidate]:
    """Run blastn -task blastn-short comparing left and right flank windows."""
    if not shutil.which("blastn"):
        return []
    candidates: list[RepeatCandidate] = []
    with tempfile.TemporaryDirectory(prefix="advena_att_") as tmpdir:
        tmp = Path(tmpdir)
        q_fa = tmp / "left.fa"
        s_fa = tmp / "right.fa"
        q_fa.write_text(f">L\n{wrap_fasta(left_seq)}\n", encoding="utf-8")
        s_fa.write_text(f">R\n{wrap_fasta(right_seq)}\n", encoding="utf-8")
        cmd = [
            "blastn",
            "-task",
            "blastn-short",
            "-word_size",
            "7",
            "-dust",
            "no",
            "-perc_identity",
            str(min_identity),
            "-strand",
            "plus",
            "-query",
            str(q_fa),
            "-subject",
            str(s_fa),
            "-outfmt",
            "6 qstart qend sstart send length pident qseq",
        ]
        res = subprocess.run(cmd, capture_output=True, text=True, check=False)
        if res.returncode != 0:
            return []
        for line in res.stdout.splitlines():
            parts = line.strip().split("\t")
            if len(parts) < 7:
                continue
            qs, qe, ss, se, alen = (int(parts[i]) for i in range(5))
            pident = float(parts[5])
            qseq = parts[6].replace("-", "")
            # Direct repeats must be on the same orientation (qs <= qe and ss <= se)
            if qs > qe or ss > se:
                continue
            if not (min_len <= alen <= max_len) or pident < min_identity:
                continue
            qseq_up = qseq.upper()
            if _is_low_complexity(qseq_up):
                continue
            candidates.append(
                RepeatCandidate(
                    attl_start=left_offset + qs - 1,
                    attl_end=left_offset + qe - 1,
                    attr_start=right_offset + ss - 1,
                    attr_end=right_offset + se - 1,
                    length=alen,
                    identity_pct=round(pident, 2),
                    seq=qseq_up,
                    method="blastn-short",
                )
            )
    return candidates


def _nearest_target(
    contig: str,
    attl_start: int,
    attl_end: int,
    attr_start: int,
    attr_end: int,
    targets_by_contig: dict[str, list[TargetFeature]],
    max_dist: int = 1000,
) -> tuple[str, int]:
    """Find closest tRNA/tmRNA/rRNA/housekeeping feature to attL or attR."""
    best_name = "-"
    best_dist = max_dist + 1
    candidates = targets_by_contig.get(contig, [])
    if not candidates:
        base = contig.split(".")[0]
        for k, v in targets_by_contig.items():
            if k.split(".")[0] == base:
                candidates = v
                break

    for feat in candidates:
        for as_, ae_ in ((attl_start, attl_end), (attr_start, attr_end)):
            if min(ae_, feat.end) >= max(as_, feat.start):
                dist = 0
            else:
                dist = min(abs(as_ - feat.end), abs(feat.start - ae_))
            if dist < best_dist:
                best_dist = dist
                best_name = feat.name
    if best_dist > max_dist:
        return ("-", -1)
    return (best_name, best_dist)


def detect_att_sites(
    sample: str,
    detected_me_tsv: Path,
    genbank_path: Path,
    gff_path: Path | None = None,
    flank_bp: int = 5000,
    min_len: int = 8,
    max_len: int = 65,
    min_identity: float = 85.0,
    use_blastn: bool = True,
) -> list[dict[str, str]]:
    """Detect attL/attR direct repeats and integration hotspots for each boundary-integrase ME."""
    if not detected_me_tsv.exists() or detected_me_tsv.stat().st_size == 0:
        return []

    records = parse_genbank(genbank_path)
    records_by_id = {r.record_id: r for r in records}
    gff_targets = _parse_targets_from_gff(gff_path)
    gbk_targets = _parse_targets_from_gbk(records_by_id)
    targets_by_contig = _merge_targets(gff_targets, gbk_targets)

    results: list[dict[str, str]] = []
    with open_text(detected_me_tsv, "rt") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            me_id = (row.get("ICE_IME_id") or "").strip()
            contig = (row.get("Genome_accession") or "").strip()
            if not me_id or me_id == "-":
                continue
            int_up = (row.get("Integrase_upstream") or "-").strip()
            int_dn = (row.get("Integrase_downstream") or "-").strip()
            has_up = bool(int_up and int_up != "-")
            has_dn = bool(int_dn and int_dn != "-")
            if not (has_up or has_dn):
                continue

            integrase_side = (
                "both" if (has_up and has_dn) else ("upstream" if has_up else "downstream")
            )
            try:
                me_start = int(float(row["Start_of_most_upstream_SP"]))
                me_end = int(float(row["Stop_of_most_downstream_SP"]))
            except (KeyError, ValueError):
                continue
            if me_start > me_end:
                me_start, me_end = me_end, me_start

            rec = records_by_id.get(contig)
            if rec is None:
                base = contig.split(".")[0]
                rec = next((r for r in records if r.accession == base or r.locus == base), None)
            if rec is None or not rec.sequence:
                continue

            seq_len = len(rec.sequence)
            # Extract L window [me_start - flank_bp, min(me_start + 1500, me_end)]
            # and R window [max(me_start, me_end - 1500), me_end + flank_bp]
            l_start = max(1, me_start - flank_bp)
            l_end = min(me_end, me_start + 1500)
            r_start = max(me_start, me_end - 1500)
            r_end = min(seq_len, me_end + flank_bp)

            # Ensure L and R windows do not overlap each other
            if l_end >= r_start:
                midpoint = (me_start + me_end) // 2
                l_end = max(l_start, midpoint)
                r_start = min(r_end, midpoint + 1)

            left_seq = rec.sequence[l_start - 1 : l_end]
            right_seq = rec.sequence[r_start - 1 : r_end]

            candidates = (
                find_direct_repeats_blastn(
                    left_seq, l_start, right_seq, r_start, min_len, max_len, min_identity
                )
                if use_blastn
                else []
            )
            if not candidates:
                candidates = find_direct_repeats_kmer(
                    left_seq, l_start, right_seq, r_start, min_len, max_len, min_identity
                )

            if not candidates:
                results.append(
                    {
                        "sample": sample,
                        "me_id": me_id,
                        "contig": rec.record_id,
                        "integrase_side": integrase_side,
                        "attl_start": "-",
                        "attl_end": "-",
                        "attr_start": "-",
                        "attr_end": "-",
                        "att_repeat_length": "0",
                        "att_identity_pct": "0.0",
                        "att_repeat_seq": "-",
                        "att_target_feature": "-",
                        "att_target_distance_bp": "-",
                        "detection_method": "none",
                    }
                )
                continue

            # Rank candidates: prioritize repeats overlapping a tRNA/housekeeping target (dist==0),
            # then adjacent targets, then longer repeat length, then higher % identity, then closer.
            scored: list[tuple[int, int, int, float, int, RepeatCandidate, str, int]] = []
            for cand in candidates:
                t_name, t_dist = _nearest_target(
                    rec.record_id,
                    cand.attl_start,
                    cand.attl_end,
                    cand.attr_start,
                    cand.attr_end,
                    targets_by_contig,
                )
                has_target = 1 if t_name != "-" else 0
                direct_overlap = 1 if (has_target and t_dist == 0) else 0
                neg_dist = -t_dist if t_dist >= 0 else -10_000
                scored.append(
                    (
                        direct_overlap,
                        has_target,
                        cand.length,
                        cand.identity_pct,
                        neg_dist,
                        cand,
                        t_name,
                        t_dist,
                    )
                )

            scored.sort(key=lambda x: (x[0], x[1], x[2], x[3], x[4]), reverse=True)
            _ov, _has_t, _len, _id, _nd, best, t_name, t_dist = scored[0]
            results.append(
                {
                    "sample": sample,
                    "me_id": me_id,
                    "contig": rec.record_id,
                    "integrase_side": integrase_side,
                    "attl_start": str(best.attl_start),
                    "attl_end": str(best.attl_end),
                    "attr_start": str(best.attr_start),
                    "attr_end": str(best.attr_end),
                    "att_repeat_length": str(best.length),
                    "att_identity_pct": f"{best.identity_pct:.2f}",
                    "att_repeat_seq": best.seq,
                    "att_target_feature": t_name,
                    "att_target_distance_bp": str(t_dist) if t_dist >= 0 else "-",
                    "detection_method": best.method,
                }
            )
    return results


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--detected-me", required=True, type=Path)
    parser.add_argument("--genbank", required=True, type=Path)
    parser.add_argument("--gff", type=Path, default=None)
    parser.add_argument("--flank-bp", type=int, default=5000)
    parser.add_argument("--min-len", type=int, default=8)
    parser.add_argument("--max-len", type=int, default=65)
    parser.add_argument("--min-identity", type=float, default=85.0)
    parser.add_argument(
        "--no-blast", action="store_true", help="Use pure-Python k-mer scanner only"
    )
    parser.add_argument("-o", "--output", required=True, type=Path)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    rows = detect_att_sites(
        sample=args.sample,
        detected_me_tsv=args.detected_me,
        genbank_path=args.genbank,
        gff_path=args.gff,
        flank_bp=args.flank_bp,
        min_len=args.min_len,
        max_len=args.max_len,
        min_identity=args.min_identity,
        use_blastn=not args.no_blast,
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with open(args.output, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=OUTPUT_COLUMNS, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    print(f"Evaluated attL/attR sites for {len(rows)} element(s) in {args.sample}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
