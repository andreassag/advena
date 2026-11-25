#!/usr/bin/env python3
"""Cluster homologous ICEs/IMEs across a cohort using skani ANI/AF and MMseqs2 protein Jaccard."""

from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

from gbk_utils import open_text

OUTPUT_COLUMNS = [
    "me_cluster_id",
    "cluster_size",
    "element_key",
    "sample",
    "me_id",
    "me_category",
    "length_bp",
    "n_cds",
    "clustering_evidence",
]


class _UnionFind:
    def __init__(self, items: list[str]) -> None:
        self.parent = {x: x for x in items}

    def find(self, x: str) -> str:
        while self.parent[x] != x:
            self.parent[x] = self.parent[self.parent[x]]
            x = self.parent[x]
        return x

    def union(self, a: str, b: str) -> None:
        ra, rb = self.find(a), self.find(b)
        if ra != rb:
            if ra < rb:
                self.parent[rb] = ra
            else:
                self.parent[ra] = rb


def _normalize_element_key(raw: str) -> str:
    """Normalize filenames or FASTA headers to '<sample>__<ME_ID>'."""
    name = Path(raw.strip()).name
    for suffix in (".fna.gz", ".fa.gz", ".fasta.gz", ".gbk.gz", ".fna", ".fa", ".fasta", ".gbk"):
        if name.endswith(suffix):
            name = name[: -len(suffix)]
            break
    return name


def _load_elements(
    bed_paths: list[Path], manifest_paths: list[Path]
) -> tuple[dict[str, dict], dict[str, str]]:
    """Load element metadata and locus_tag -> element_key mapping."""
    elements: dict[str, dict] = {}
    locus_to_elem: dict[str, str] = {}

    for bed in bed_paths:
        if not bed.exists() or bed.stat().st_size == 0:
            continue
        sample = bed.name.removesuffix("_ME_regions.bed")
        with open_text(bed, "rt") as handle:
            for line in handle:
                if not line.strip() or line.startswith("#"):
                    continue
                parts = line.rstrip("\r\n").split("\t")
                if len(parts) < 4:
                    continue
                start, end = int(parts[1]), int(parts[2])
                me_id = parts[3].strip()
                cat = parts[6].strip() if len(parts) >= 7 else "Unknown"
                key = f"{sample}__{me_id}"
                elements[key] = {
                    "sample": sample,
                    "me_id": me_id,
                    "me_category": cat,
                    "length_bp": max(0, end - start),
                    "n_cds": 0,
                }

    for manifest in manifest_paths:
        if not manifest.exists() or manifest.stat().st_size == 0:
            continue
        with open_text(manifest, "rt") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                sample = (row.get("sample") or "").strip()
                me_id = (row.get("me_id") or "").strip()
                lt = (row.get("locus_tag") or "").strip()
                if not sample or not me_id or me_id == "GENOME_BG":
                    continue
                key = f"{sample}__{me_id}"
                info = elements.setdefault(
                    key,
                    {
                        "sample": sample,
                        "me_id": me_id,
                        "me_category": row.get("me_category", "Unknown"),
                        "length_bp": 0,
                        "n_cds": 0,
                    },
                )
                info["n_cds"] += 1
                if lt:
                    locus_to_elem[lt] = key
                    locus_to_elem[f"{sample}|{lt}"] = key

    return elements, locus_to_elem


def _parse_skani_edges(
    skani_path: Path | None, ani_threshold: float, af_threshold: float
) -> list[tuple[str, str]]:
    """Parse skani triangle --sparse or skani dist TSV into qualifying edges."""
    if skani_path is None or not skani_path.exists() or skani_path.stat().st_size == 0:
        return []
    edges: list[tuple[str, str]] = []
    with open_text(skani_path, "rt") as handle:
        first = handle.readline().strip()
        if not first:
            return []
        if "Ref_file" in first or "ANI" in first:
            handle.seek(0)
            for row in csv.DictReader(handle, delimiter="\t"):
                ref = _normalize_element_key(row.get("Ref_file") or row.get("Ref_name") or "")
                qry = _normalize_element_key(row.get("Query_file") or row.get("Query_name") or "")
                if not ref or not qry or ref == qry:
                    continue
                try:
                    ani = float(row.get("ANI") or 0.0)
                    af_ref = float(row.get("Align_fraction_ref") or 0.0)
                    af_qry = float(row.get("Align_fraction_query") or 0.0)
                except ValueError:
                    continue
                if ani >= ani_threshold and min(af_ref, af_qry) >= af_threshold:
                    edges.append((ref, qry))
        else:
            # Lower-triangular PHYLIP matrix from skani triangle
            lines = [ln.strip() for ln in handle if ln.strip()]
            names: list[str] = []
            for line in lines:
                parts = line.split("\t")
                name = _normalize_element_key(parts[0])
                names.append(name)
                for j, val in enumerate(parts[1:]):
                    try:
                        ani = float(val)
                    except ValueError:
                        continue
                    if ani >= ani_threshold and j < len(names) - 1:
                        edges.append((names[j], name))
    return edges


def _parse_mmseqs_jaccard_edges(
    mmseqs_path: Path | None,
    locus_to_elem: dict[str, str],
    elements: dict[str, dict],
    min_jaccard: float = 0.50,
    short_len_bp: int = 5000,
) -> list[tuple[str, str]]:
    """Build protein-cluster Jaccard edges for short IMEs (<short_len_bp) or missing-length MEs."""
    if mmseqs_path is None or not mmseqs_path.exists() or mmseqs_path.stat().st_size == 0:
        return []

    elem_clusters: dict[str, set[str]] = {}
    with open_text(mmseqs_path, "rt") as handle:
        for line in handle:
            parts = line.rstrip("\r\n").split("\t")
            if len(parts) < 2:
                continue
            rep, member = parts[0].strip(), parts[1].strip()
            elem = locus_to_elem.get(member)
            if elem:
                elem_clusters.setdefault(elem, set()).add(rep)

    short_elems = sorted(
        k
        for k, info in elements.items()
        if (info["length_bp"] == 0 or info["length_bp"] < short_len_bp) and k in elem_clusters
    )
    edges: list[tuple[str, str]] = []
    for i in range(len(short_elems)):
        for j in range(i + 1, len(short_elems)):
            a, b = short_elems[i], short_elems[j]
            sa, sb = elem_clusters[a], elem_clusters[b]
            union_size = len(sa | sb)
            if union_size == 0:
                continue
            jaccard = len(sa & sb) / union_size
            if jaccard >= min_jaccard:
                edges.append((a, b))
    return edges


def cluster_elements(
    bed_paths: list[Path],
    manifest_paths: list[Path],
    skani_path: Path | None = None,
    mmseqs_path: Path | None = None,
    ani_threshold: float = 95.0,
    af_threshold: float = 70.0,
    protein_jaccard: float = 0.50,
) -> tuple[list[dict[str, str]], dict[str, list[str]]]:
    """Assign deterministic CLUSTER_001..N IDs via Union-Find."""
    elements, locus_to_elem = _load_elements(bed_paths, manifest_paths)

    skani_edges = _parse_skani_edges(skani_path, ani_threshold, af_threshold)
    for a, b in skani_edges:
        for k in (a, b):
            if k not in elements:
                sample, me_id = k.split("__", 1) if "__" in k else (k, "ME_1")
                elements[k] = {
                    "sample": sample,
                    "me_id": me_id,
                    "me_category": "Unknown",
                    "length_bp": 0,
                    "n_cds": 0,
                }

    mmseqs_edges = _parse_mmseqs_jaccard_edges(
        mmseqs_path, locus_to_elem, elements, min_jaccard=protein_jaccard
    )

    if not elements:
        return ([], {})

    uf = _UnionFind(sorted(elements))
    evidence: dict[str, set[str]] = {k: set() for k in elements}
    for a, b in skani_edges:
        if a in elements and b in elements:
            uf.union(a, b)
            evidence[a].add("skani_ani")
            evidence[b].add("skani_ani")
    for a, b in mmseqs_edges:
        if a in elements and b in elements:
            uf.union(a, b)
            evidence[a].add("mmseqs_jaccard")
            evidence[b].add("mmseqs_jaccard")

    raw_clusters: dict[str, list[str]] = {}
    for k in sorted(elements):
        raw_clusters.setdefault(uf.find(k), []).append(k)

    ordered_clusters = sorted(
        raw_clusters.values(), key=lambda members: (-len(members), members[0])
    )

    rows: list[dict[str, str]] = []
    cluster_members: dict[str, list[str]] = {}
    for idx, members in enumerate(ordered_clusters, start=1):
        cid = f"CLUSTER_{idx:03d}"
        cluster_members[cid] = members
        for key in members:
            info = elements[key]
            ev = ";".join(sorted(evidence[key])) if evidence[key] else "singleton"
            rows.append(
                {
                    "me_cluster_id": cid,
                    "cluster_size": str(len(members)),
                    "element_key": key,
                    "sample": info["sample"],
                    "me_id": info["me_id"],
                    "me_category": info["me_category"],
                    "length_bp": str(info["length_bp"]),
                    "n_cds": str(info["n_cds"]),
                    "clustering_evidence": ev,
                }
            )
    return rows, cluster_members


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--beds", nargs="*", type=Path, default=[])
    parser.add_argument("--manifests", nargs="*", type=Path, default=[])
    parser.add_argument("--skani", type=Path, default=None)
    parser.add_argument("--mmseqs", type=Path, default=None)
    parser.add_argument("--ani-threshold", type=float, default=95.0)
    parser.add_argument("--af-threshold", type=float, default=70.0)
    parser.add_argument("--protein-jaccard", type=float, default=0.50)
    parser.add_argument("--out-tsv", required=True, type=Path)
    parser.add_argument("--out-json", required=True, type=Path)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    rows, members = cluster_elements(
        bed_paths=args.beds,
        manifest_paths=args.manifests,
        skani_path=args.skani,
        mmseqs_path=args.mmseqs,
        ani_threshold=args.ani_threshold,
        af_threshold=args.af_threshold,
        protein_jaccard=args.protein_jaccard,
    )
    args.out_tsv.parent.mkdir(parents=True, exist_ok=True)
    with open(args.out_tsv, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=OUTPUT_COLUMNS, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    args.out_json.parent.mkdir(parents=True, exist_ok=True)
    args.out_json.write_text(json.dumps(members, indent=2) + "\n", encoding="utf-8")
    print(f"Clustered {len(rows)} element(s) into {len(members)} cluster(s)", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
