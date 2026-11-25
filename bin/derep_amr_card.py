#!/usr/bin/env python3
"""Native locus-tag and coordinate-overlap dereplication of AMRFinderPlus (--plus) and CARD RGI."""

from __future__ import annotations

import argparse
import csv
import re
import sys
from dataclasses import dataclass
from pathlib import Path

from gbk_utils import open_text


@dataclass
class HitRecord:
    """Normalized AMR/stress/biocide/virulence hit from AMRFinderPlus or CARD RGI."""

    source: str  # "amrfinderplus" | "card_rgi"
    sample: str
    me_id: str
    contig: str
    locus_tag: str
    start: int
    end: int
    strand: str
    symbol: str
    name: str
    element_type: str
    element_subtype: str
    drug_classes: list[str]
    method_or_cutoff: str
    pct_identity: float
    pct_coverage: float
    tier: int
    aro_accession: str = ""
    aro_name: str = ""
    resistance_mechanism: str = ""
    accession: str = ""


class _UnionFind:
    def __init__(self, n: int) -> None:
        self.parent = list(range(n))

    def find(self, x: int) -> int:
        while self.parent[x] != x:
            self.parent[x] = self.parent[self.parent[x]]
            x = self.parent[x]
        return x

    def union(self, a: int, b: int) -> None:
        ra, rb = self.find(a), self.find(b)
        if ra != rb:
            if ra < rb:
                self.parent[rb] = ra
            else:
                self.parent[ra] = rb


def _safe_float(val: str | None, default: float = 0.0) -> float:
    if val is None:
        return default
    cleaned = val.strip().rstrip("%")
    if not cleaned or cleaned in {"-", "NA", "N/A"}:
        return default
    try:
        return float(cleaned)
    except ValueError:
        return default


def _safe_int(val: str | None, default: int = 0) -> int:
    if val is None:
        return default
    cleaned = val.strip()
    if not cleaned or cleaned in {"-", "NA", "N/A"}:
        return default
    try:
        return int(float(cleaned))
    except ValueError:
        return default


def _normalize_strand(val: str | None) -> str:
    if not val:
        return "."
    s = val.strip()
    if s in {"+", "1", "forward", "fwd"}:
        return "+"
    if s in {"-", "-1", "reverse", "rev"}:
        return "-"
    return "."


def _normalize_drug_tokens(raw_values: list[str]) -> list[str]:
    tokens: list[str] = []
    seen: set[str] = set()
    for raw in raw_values:
        if not raw or raw.strip() in {"-", "NA", "N/A"}:
            continue
        for part in re.split(r"[;/,]", raw):
            norm = re.sub(r"\s+", " ", part.strip()).lower()
            norm = re.sub(r"\s+antibiotic$", "", norm).strip()
            if norm and norm not in seen:
                seen.add(norm)
                tokens.append(norm)
    return sorted(tokens)


def _amrfinder_tier(method: str) -> int:
    m = method.upper().strip()
    if any(m.startswith(tok) for tok in ("ALLELE", "EXACT")):
        return 3
    if "PARTIAL" in m or m.startswith("HMM") or "INTERNAL_STOP" in m:
        return 1
    if any(m.startswith(tok) for tok in ("BLAST", "POINT")):
        return 2
    return 1


def _rgi_tier(cutoff: str) -> int:
    c = cutoff.strip().lower()
    if c == "perfect":
        return 3
    if c == "strict":
        return 2
    return 1


def reciprocal_overlap(start_a: int, end_a: int, start_b: int, end_b: int) -> float:
    """Compute reciprocal coordinate overlap fraction for 1-based closed intervals."""
    if start_a <= 0 or end_a < start_a or start_b <= 0 or end_b < start_b:
        return 0.0
    ov = max(0, min(end_a, end_b) - max(start_a, start_b) + 1)
    if ov == 0:
        return 0.0
    len_a = end_a - start_a + 1
    len_b = end_b - start_b + 1
    return min(ov / len_a, ov / len_b)


def compatible_strand(strand_a: str, strand_b: str) -> bool:
    """Return True if strands match or either strand is unspecified."""
    if strand_a in {"", "."} or strand_b in {"", "."}:
        return True
    return strand_a == strand_b


def load_manifest(manifest_path: Path | None) -> dict[str, dict[str, str]]:
    """Load *_ME_cds_manifest.tsv indexed by locus_tag."""
    by_locus: dict[str, dict[str, str]] = {}
    if manifest_path is None or not manifest_path.exists() or manifest_path.stat().st_size == 0:
        return by_locus
    with open_text(manifest_path, "rt") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            lt = (row.get("locus_tag") or "").strip()
            if lt:
                by_locus[lt] = row
    return by_locus


def load_me_bed(bed_path: Path | None) -> dict[str, tuple[str, int, int]]:
    """Load *_ME_regions.bed indexed by me_id -> (contig, 1-based slice_start, slice_end)."""
    by_me: dict[str, tuple[str, int, int]] = {}
    if bed_path is None or not bed_path.exists() or bed_path.stat().st_size == 0:
        return by_me
    with open_text(bed_path, "rt") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\r\n").split("\t")
            if len(parts) < 4:
                continue
            try:
                s0, e1 = int(parts[1]), int(parts[2])
            except ValueError:
                continue
            contig = parts[0].strip()
            me_id = parts[3].strip()
            if contig and me_id:
                by_me[me_id] = (contig, s0 + 1, e1)
    return by_me


def parse_amrfinder(
    path: Path | None,
    sample: str,
    manifest: dict[str, dict[str, str]],
    me_bed: dict[str, tuple[str, int, int]] | None = None,
) -> list[HitRecord]:
    """Parse AMRFinderPlus (--plus) TSV, keeping AMR, STRESS, METAL, BIOCIDE, and VIRULENCE."""
    if path is None or not path.exists() or path.stat().st_size == 0:
        return []
    me_bed = me_bed or {}
    hits: list[HitRecord] = []
    with open_text(path, "rt") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            prot_id = (
                row.get("Protein id") or row.get("Protein identifier") or row.get("Name") or ""
            ).strip()
            locus_tag = prot_id.split()[0] if prot_id and prot_id != "NA" else ""
            m_info = manifest.get(locus_tag, {})

            raw_contig = (row.get("Contig id") or "").strip()
            contig = raw_contig
            raw_start = _safe_int(row.get("Start"))
            raw_end = _safe_int(row.get("Stop"))
            if raw_start > raw_end and raw_end > 0:
                raw_start, raw_end = raw_end, raw_start

            me_id = m_info.get("me_id", "UNASSIGNED")
            lifted_contig = ""
            lifted_start = 0
            lifted_end = 0
            if not m_info and raw_contig and "__" in raw_contig:
                prefix_part, cand_me = raw_contig.rsplit("__", 1)
                if prefix_part == sample or cand_me in me_bed:
                    me_id = cand_me
                    if cand_me in me_bed:
                        orig_contig, slice_s, _slice_e = me_bed[cand_me]
                        lifted_contig = orig_contig
                        if raw_start > 0 and raw_end > 0:
                            lifted_start = slice_s + raw_start - 1
                            lifted_end = slice_s + raw_end - 1

            if m_info.get("contig"):
                contig = m_info["contig"]
            elif lifted_contig:
                contig = lifted_contig
            elif not contig or contig == "NA":
                contig = "UNKNOWN"

            start = _safe_int(m_info.get("orig_start")) or lifted_start or raw_start
            end = _safe_int(m_info.get("orig_end")) or lifted_end or raw_end
            if start > end and end > 0:
                start, end = end, start
            strand = _normalize_strand(m_info.get("strand") or row.get("Strand"))

            symbol = (row.get("Element symbol") or row.get("Gene symbol") or "").strip()
            name = (row.get("Element name") or row.get("Sequence name") or "").strip()
            elem_type = (row.get("Type") or row.get("Element type") or "AMR").strip().upper()
            elem_subtype = (
                (row.get("Subtype") or row.get("Element subtype") or elem_type).strip().upper()
            )
            if not symbol and not name:
                continue

            drug_cls = _normalize_drug_tokens([row.get("Class") or "", row.get("Subclass") or ""])
            method = (row.get("Method") or "").strip()
            pct_id = _safe_float(
                row.get("% Identity to reference") or row.get("% Identity to reference sequence")
            )
            pct_cov = _safe_float(
                row.get("% Coverage of reference") or row.get("% Coverage of reference sequence")
            )
            acc = (
                row.get("Closest reference accession")
                or row.get("Accession of closest sequence")
                or ""
            ).strip()

            hits.append(
                HitRecord(
                    source="amrfinderplus",
                    sample=sample,
                    me_id=me_id,
                    contig=contig,
                    locus_tag=locus_tag,
                    start=start,
                    end=end,
                    strand=strand,
                    symbol=symbol or name,
                    name=name or symbol,
                    element_type=elem_type,
                    element_subtype=elem_subtype,
                    drug_classes=drug_cls,
                    method_or_cutoff=method or "UNKNOWN",
                    pct_identity=pct_id,
                    pct_coverage=pct_cov,
                    tier=_amrfinder_tier(method),
                    accession=acc,
                )
            )
    return hits


def parse_card_rgi(
    path: Path | None,
    sample: str,
    manifest: dict[str, dict[str, str]],
    include_loose: bool = False,
) -> list[HitRecord]:
    """Parse CARD RGI main TSV output (Perfect + Strict by default)."""
    if path is None or not path.exists() or path.stat().st_size == 0:
        return []
    allowed = {"perfect", "strict"}
    if include_loose:
        allowed.add("loose")

    hits: list[HitRecord] = []
    with open_text(path, "rt") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            cutoff = (row.get("Cut_Off") or "").strip()
            if cutoff.lower() not in allowed:
                continue

            orf_id = (row.get("ORF_ID") or "").strip()
            raw_token = orf_id.split()[0] if orf_id else ""
            locus_tag = raw_token if raw_token in manifest else ""
            if not locus_tag and raw_token:
                # Try stripping trailing _1 added by prodigal/RGI contig mode
                stripped = re.sub(r"_\d+$", "", raw_token)
                locus_tag = stripped if stripped in manifest else raw_token

            m_info = manifest.get(locus_tag, {})
            contig = m_info.get("contig") or (row.get("Contig") or "").strip() or "UNKNOWN"
            start = _safe_int(m_info.get("orig_start")) or _safe_int(row.get("Start"))
            end = _safe_int(m_info.get("orig_end")) or _safe_int(row.get("Stop"))
            if start > end and end > 0:
                start, end = end, start
            strand = _normalize_strand(m_info.get("strand") or row.get("Orientation"))
            me_id = m_info.get("me_id", "UNASSIGNED")

            aro_num = (row.get("ARO") or "").strip()
            if aro_num and not aro_num.upper().startswith("ARO:"):
                aro_acc = f"ARO:{aro_num}"
            else:
                aro_acc = aro_num
            aro_name = (row.get("Best_Hit_ARO") or "").strip()
            mech = (row.get("Resistance Mechanism") or "").strip()
            drug_cls = _normalize_drug_tokens([row.get("Drug Class") or ""])
            pct_id = _safe_float(row.get("Best_Identities"))
            pct_cov = _safe_float(row.get("Percentage Length of Reference Sequence"))

            hits.append(
                HitRecord(
                    source="card_rgi",
                    sample=sample,
                    me_id=me_id,
                    contig=contig,
                    locus_tag=locus_tag,
                    start=start,
                    end=end,
                    strand=strand,
                    symbol=aro_name or aro_acc,
                    name=aro_name or aro_acc,
                    element_type="AMR",
                    element_subtype="AMR",
                    drug_classes=drug_cls,
                    method_or_cutoff=cutoff,
                    pct_identity=pct_id,
                    pct_coverage=pct_cov,
                    tier=_rgi_tier(cutoff),
                    aro_accession=aro_acc,
                    aro_name=aro_name,
                    resistance_mechanism=mech,
                    accession=aro_acc,
                )
            )
    return hits


def _rank_key(hit: HitRecord) -> tuple[int, float, int]:
    source_priority = 2 if hit.source == "amrfinderplus" else 1
    return (hit.tier, round(hit.pct_identity * hit.pct_coverage, 2), source_priority)


def dereplicate_hits(hits: list[HitRecord], min_overlap: float = 0.80) -> list[dict[str, str]]:
    """Cluster hits by locus_tag or >=min_overlap reciprocal overlap and merge fields."""
    if not hits:
        return []

    by_contig: dict[tuple[str, str], list[HitRecord]] = {}
    for h in hits:
        by_contig.setdefault((h.sample, h.contig), []).append(h)

    merged_rows: list[dict[str, str]] = []

    for (sample, contig), group in sorted(by_contig.items()):
        uf = _UnionFind(len(group))
        for i in range(len(group)):
            for j in range(i + 1, len(group)):
                a, b = group[i], group[j]
                same_locus = bool(a.locus_tag and b.locus_tag and a.locus_tag == b.locus_tag)
                coord_match = (
                    compatible_strand(a.strand, b.strand)
                    and reciprocal_overlap(a.start, a.end, b.start, b.end) >= min_overlap
                )
                if same_locus or coord_match:
                    uf.union(i, j)

        clusters: dict[int, list[HitRecord]] = {}
        for idx, hit in enumerate(group):
            clusters.setdefault(uf.find(idx), []).append(hit)

        for cluster in clusters.values():
            ranked = sorted(cluster, key=_rank_key, reverse=True)
            rep = ranked[0]
            amr_hits = [h for h in ranked if h.source == "amrfinderplus"]
            card_hits = [h for h in ranked if h.source == "card_rgi"]
            best_amr = amr_hits[0] if amr_hits else None
            best_card = card_hits[0] if card_hits else None

            if best_amr and best_card:
                detected_by = "both"
            elif best_amr:
                detected_by = "amrfinderplus"
            else:
                detected_by = "card_rgi"

            # Canonical nomenclature prefers AMRFinderPlus curated bacterial symbols when available
            canonical_symbol = best_amr.symbol if best_amr else rep.symbol
            element_type = best_amr.element_type if best_amr else rep.element_type
            element_subtype = best_amr.element_subtype if best_amr else rep.element_subtype

            all_classes = _normalize_drug_tokens([cls for h in cluster for cls in h.drug_classes])
            locus_tag = next((h.locus_tag for h in ranked if h.locus_tag), "-")
            me_id = next((h.me_id for h in ranked if h.me_id != "UNASSIGNED"), "UNASSIGNED")
            starts = [h.start for h in cluster if h.start > 0]
            ends = [h.end for h in cluster if h.end > 0]
            strand = next((h.strand for h in ranked if h.strand in {"+", "-"}), ".")

            merged_rows.append(
                {
                    "sample": sample,
                    "me_id": me_id,
                    "contig": contig,
                    "locus_tag": locus_tag,
                    "start": str(min(starts)) if starts else "0",
                    "end": str(max(ends)) if ends else "0",
                    "strand": strand,
                    "canonical_symbol": canonical_symbol,
                    "element_type": element_type,
                    "element_subtype": element_subtype,
                    "card_aro_accession": best_card.aro_accession if best_card else "-",
                    "card_aro_name": best_card.aro_name if best_card else "-",
                    "resistance_mechanism": (
                        best_card.resistance_mechanism
                        if (best_card and best_card.resistance_mechanism)
                        else "-"
                    ),
                    "consensus_drug_class": ";".join(all_classes) if all_classes else "-",
                    "detected_by": detected_by,
                    "representative_source": rep.source,
                    "representative_tier": str(rep.tier),
                    "amrfinder_symbol": best_amr.symbol if best_amr else "-",
                    "amrfinder_method": best_amr.method_or_cutoff if best_amr else "-",
                    "amrfinder_identity": f"{best_amr.pct_identity:.2f}" if best_amr else "-",
                    "amrfinder_coverage": f"{best_amr.pct_coverage:.2f}" if best_amr else "-",
                    "card_cut_off": best_card.method_or_cutoff if best_card else "-",
                    "card_identity": f"{best_card.pct_identity:.2f}" if best_card else "-",
                    "card_coverage": f"{best_card.pct_coverage:.2f}" if best_card else "-",
                }
            )

    merged_rows.sort(
        key=lambda r: (
            r["sample"],
            r["contig"],
            int(r["start"]),
            r["locus_tag"],
            r["canonical_symbol"],
        )
    )
    return merged_rows


OUTPUT_COLUMNS = [
    "sample",
    "me_id",
    "contig",
    "locus_tag",
    "start",
    "end",
    "strand",
    "canonical_symbol",
    "element_type",
    "element_subtype",
    "card_aro_accession",
    "card_aro_name",
    "resistance_mechanism",
    "consensus_drug_class",
    "detected_by",
    "representative_source",
    "representative_tier",
    "amrfinder_symbol",
    "amrfinder_method",
    "amrfinder_identity",
    "amrfinder_coverage",
    "card_cut_off",
    "card_identity",
    "card_coverage",
]


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sample", required=True, help="Sample identifier")
    parser.add_argument("--manifest", type=Path, default=None, help="*_ME_cds_manifest.tsv")
    parser.add_argument("--bed", type=Path, default=None, help="*_ME_regions.bed for slice lifting")
    parser.add_argument("--amrfinder", type=Path, default=None, help="AMRFinderPlus TSV")
    parser.add_argument("--card-rgi", type=Path, default=None, help="CARD RGI main TXT")
    parser.add_argument("--min-overlap", type=float, default=0.80, help="Reciprocal overlap cutoff")
    parser.add_argument("--include-loose", action="store_true", help="Include CARD RGI Loose hits")
    parser.add_argument("-o", "--output", required=True, type=Path, help="Dereplicated output TSV")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    manifest = load_manifest(args.manifest)
    me_bed = load_me_bed(args.bed)
    hits = parse_amrfinder(args.amrfinder, args.sample, manifest, me_bed=me_bed) + parse_card_rgi(
        args.card_rgi, args.sample, manifest, include_loose=args.include_loose
    )
    rows = dereplicate_hits(hits, min_overlap=args.min_overlap)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with open(args.output, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=OUTPUT_COLUMNS, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    print(
        f"Dereplicated {len(hits)} raw AMR/CARD hit(s) into {len(rows)} non-redundant call(s) "
        f"for {args.sample}",
        file=sys.stderr,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
