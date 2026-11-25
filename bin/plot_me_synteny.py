#!/usr/bin/env python3
"""Render SVG and interactive HTML synteny alignments for clustered ICE/IME GenBank slices."""

from __future__ import annotations

import argparse
import csv
import html
import sys
from pathlib import Path

from gbk_utils import open_text, parse_genbank

COLOR_PALETTE = {
    "ICE_BACKBONE_SP": "#1f77b4",  # Blue: Signature Proteins
    "AMR_CARD": "#d62728",  # Red: AMR / Stress / Biocide / CARD
    "DEFENSE_TA": "#9467bd",  # Purple: DefenseFinder / TA
    "BGC": "#2ca02c",  # Green: Secondary metabolite BGC
    "CARGO_CHARACTERIZED_OG": "#ff7f0e",  # Orange: Characterized eggNOG OG
    "CARGO_HYPOTHETICAL_OG": "#bcbd22",  # Olive/Grey: Hypothetical / uncharacterized
}


def _load_cargo_roles(cargo_paths: list[Path]) -> dict[tuple[str, str], tuple[str, str]]:
    """Map (sample, locus_tag) -> (functional_layer, label)."""
    roles: dict[tuple[str, str], tuple[str, str]] = {}
    for p in cargo_paths:
        if not p.exists() or p.stat().st_size == 0:
            continue
        with open_text(p, "rt") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                sample = (row.get("sample") or "").strip()
                lt = (row.get("locus_tag") or "").strip()
                if not sample or not lt:
                    continue
                if row.get("is_signature_protein", "").upper() == "TRUE":
                    layer = "ICE_BACKBONE_SP"
                    label = row.get("sp_type") or lt
                elif (row.get("amr_gene") or "-") != "-":
                    layer = "AMR_CARD"
                    label = row["amr_gene"]
                elif (row.get("defense_system") or "-") != "-":
                    layer = "DEFENSE_TA"
                    label = row["defense_system"]
                elif (row.get("bgc_type") or "-") != "-":
                    layer = "BGC"
                    label = row["bgc_type"]
                elif row.get("cargo_role") == "CARGO_CHARACTERIZED_OG":
                    layer = "CARGO_CHARACTERIZED_OG"
                    label = row.get("eggnog_bacillota_og") or lt
                else:
                    layer = "CARGO_HYPOTHETICAL_OG"
                    label = lt
                roles[(sample, lt)] = (layer, label)
    return roles


def _index_gbk_slices(gbk_paths: list[Path]) -> dict[str, Path]:
    index: dict[str, Path] = {}
    for p in gbk_paths:
        name = p.name
        for suf in (".gbk.gz", ".gb.gz", ".gbff.gz", ".gbk", ".gb", ".gbff"):
            if name.endswith(suf):
                name = name[: -len(suf)]
                break
        index[name] = p
    return index


def _render_standalone_svg(
    cluster_id: str,
    member_records: list[tuple[str, str, object]],
    cargo_roles: dict[tuple[str, str], tuple[str, str]],
) -> str:
    """Render a clean publication-ready SVG gene track alignment with legend."""
    width = 1100
    margin_left = 190
    margin_right = 40
    track_h = 80
    top_pad = 65
    plot_w = width - margin_left - margin_right
    max_len = max((len(rec.sequence) for _, _, rec in member_records), default=1000) or 1000
    height = top_pad + len(member_records) * track_h + 70

    svg: list[str] = [
        f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {width} {height}" '
        f'width="{width}" height="{height}" font-family="DejaVu Sans, Arial, sans-serif">',
        f'<rect x="0" y="0" width="{width}" height="{height}" fill="#ffffff"/>',
        f'<text x="20" y="32" font-size="18" font-weight="bold" fill="#1a1a1a">'
        f"advena Synteny Alignment — {html.escape(cluster_id)} "
        f"({len(member_records)} elements)</text>",
    ]

    for row_idx, (elem_key, sample, rec) in enumerate(member_records):
        y_mid = top_pad + row_idx * track_h + 32
        seq_len = max(1, len(rec.sequence))
        bar_w = max(20.0, plot_w * (seq_len / max_len))

        svg.append(
            f'<text x="20" y="{y_mid + 4}" font-size="12" font-weight="600" fill="#222">'
            f"{html.escape(elem_key)}</text>"
        )
        svg.append(
            f'<text x="20" y="{y_mid + 19}" font-size="10" fill="#666">{seq_len:,} bp</text>'
        )
        svg.append(
            f'<line x1="{margin_left}" y1="{y_mid}" x2="{margin_left + bar_w:.1f}" '
            f'y2="{y_mid}" stroke="#999" stroke-width="2"/>'
        )

        for feat in rec.features:
            if feat.feature_type != "CDS":
                continue
            lt = feat.get_first("locus_tag") or feat.get_first("protein_id") or ""
            notes = " ".join(feat.qualifiers.get("note", []))
            if (sample, lt) in cargo_roles:
                layer, label = cargo_roles[(sample, lt)]
            elif "ICEscreen_SP:" in notes:
                layer = "ICE_BACKBONE_SP"
                label = notes.split("ICEscreen_SP:", 1)[1].split()[0]
            else:
                layer = "CARGO_HYPOTHETICAL_OG"
                label = lt

            color = COLOR_PALETTE.get(layer, "#bcbd22")
            x1 = margin_left + (feat.start / max_len) * plot_w
            x2 = margin_left + (feat.end / max_len) * plot_w
            w = max(4.0, x2 - x1)
            head = min(8.0, w * 0.4)
            h = 10.0

            if feat.strand == "-":
                pts = (
                    f"{x1:.1f},{y_mid} {x1 + head:.1f},{y_mid - h} {x2:.1f},{y_mid - h} "
                    f"{x2:.1f},{y_mid + h} {x1 + head:.1f},{y_mid + h}"
                )
            else:
                pts = (
                    f"{x1:.1f},{y_mid - h} {x2 - head:.1f},{y_mid - h} {x2:.1f},{y_mid} "
                    f"{x2 - head:.1f},{y_mid + h} {x1:.1f},{y_mid + h}"
                )
            title = html.escape(f"{lt} ({layer}: {label}) [{feat.start}..{feat.end} {feat.strand}]")
            svg.append(
                f'<polygon points="{pts}" fill="{color}" stroke="#222" stroke-width="0.6">'
                f"<title>{title}</title></polygon>"
            )
            if layer in {"ICE_BACKBONE_SP", "AMR_CARD", "DEFENSE_TA", "BGC"} and w >= 18:
                svg.append(
                    f'<text x="{(x1 + x2) / 2:.1f}" y="{y_mid - 14}" font-size="8.5" '
                    f'text-anchor="middle" fill="#111">{html.escape(label[:14])}</text>'
                )

    # Legend
    ly = height - 34
    lx = 20
    for layer, color in COLOR_PALETTE.items():
        svg.append(f'<rect x="{lx}" y="{ly}" width="12" height="12" fill="{color}" stroke="#333"/>')
        svg.append(
            f'<text x="{lx + 16}" y="{ly + 10}" font-size="10" fill="#222">'
            f"{html.escape(layer)}</text>"
        )
        lx += 175

    svg.append("</svg>\n")
    return "\n".join(svg)


def _render_with_pygenomeviz(
    cluster_id: str,
    gbk_files: list[tuple[str, str, Path]],
    cargo_roles: dict[tuple[str, str], tuple[str, str]],
    out_svg: Path,
    out_html: Path,
) -> bool:
    """Use pyGenomeViz when installed; return False to trigger standalone SVG fallback."""
    try:
        from pygenomeviz import GenomeViz  # type: ignore[import-not-found]
        from pygenomeviz.parser import Genbank  # type: ignore[import-not-found]
    except ImportError:
        return False

    try:
        gv = GenomeViz(fig_track_height=0.7, feature_track_ratio=0.35)
        gv.set_scale_bar()
        for elem_key, sample, gbk_path in gbk_files:
            gbk = Genbank(gbk_path)
            track = gv.add_feature_track(elem_key, gbk.get_seqid2size())
            for seqid, features in gbk.get_seqid2features("CDS").items():
                segment = track.get_segment(seqid)
                for feat in features:
                    lt = (feat.qualifiers.get("locus_tag") or [""])[0] or (
                        feat.qualifiers.get("protein_id") or [""]
                    )[0]
                    notes = " ".join(feat.qualifiers.get("note", []))
                    if (sample, lt) in cargo_roles:
                        layer, label = cargo_roles[(sample, lt)]
                    elif "ICEscreen_SP:" in notes:
                        layer = "ICE_BACKBONE_SP"
                        label = notes.split("ICEscreen_SP:", 1)[1].split()[0]
                    else:
                        layer = "CARGO_HYPOTHETICAL_OG"
                        label = lt
                    color = COLOR_PALETTE.get(layer, "#bcbd22")
                    show_label = (
                        label
                        if layer in {"ICE_BACKBONE_SP", "AMR_CARD", "DEFENSE_TA", "BGC"}
                        else ""
                    )
                    segment.add_feature(
                        int(feat.location.start),
                        int(feat.location.end),
                        int(feat.location.strand or 1),
                        label=show_label,
                        facecolor=color,
                        plotstyle="bigarrow",
                    )
        gv.savefig(out_svg)
        gv.savefig_html(out_html)
        return True
    except Exception as exc:  # noqa: BLE001
        print(
            f"WARNING: pyGenomeViz rendering failed for {cluster_id} ({exc}); "
            "falling back to standalone SVG.",
            file=sys.stderr,
        )
        return False


def plot_clusters(
    clusters_tsv: Path,
    gbk_slices: list[Path],
    cargo_profiles: list[Path],
    out_dir: Path,
    max_elements: int = 15,
) -> int:
    """Generate SVG and HTML synteny plots for each multi-member cluster."""
    out_dir.mkdir(parents=True, exist_ok=True)
    if not clusters_tsv.exists() or clusters_tsv.stat().st_size == 0:
        return 0

    gbk_index = _index_gbk_slices(gbk_slices)
    cargo_roles = _load_cargo_roles(cargo_profiles)

    clusters: dict[str, list[tuple[str, str]]] = {}
    with open_text(clusters_tsv, "rt") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            cid = (row.get("me_cluster_id") or "").strip()
            key = (row.get("element_key") or "").strip()
            sample = (row.get("sample") or "").strip()
            if cid and key:
                clusters.setdefault(cid, []).append((key, sample))

    plotted = 0
    for cid, members in sorted(clusters.items()):
        if len(members) < 2:
            continue
        selected = members[:max_elements]
        resolved = [(k, s, gbk_index[k]) for k, s in selected if k in gbk_index]
        if len(resolved) < 2:
            continue

        out_svg = out_dir / f"{cid}_synteny.svg"
        out_html = out_dir / f"{cid}_synteny.html"

        used_pgv = _render_with_pygenomeviz(cid, resolved, cargo_roles, out_svg, out_html)
        if not used_pgv:
            parsed_members = []
            for key, sample, gbk_path in resolved:
                recs = parse_genbank(gbk_path)
                if recs:
                    parsed_members.append((key, sample, recs[0]))
            if len(parsed_members) < 2:
                continue
            svg_content = _render_standalone_svg(cid, parsed_members, cargo_roles)
            out_svg.write_text(svg_content, encoding="utf-8")
            html_content = (
                "<!DOCTYPE html>\n<html lang='en'>\n<head>\n<meta charset='utf-8'>\n"
                f"<title>advena Synteny — {html.escape(cid)}</title>\n"
                "<style>body{font-family:system-ui,sans-serif;margin:24px;background:#fafafa}"
                ".card{background:#fff;border:1px solid #ddd;border-radius:8px;padding:16px;"
                "overflow-x:auto}</style>\n</head>\n<body>\n"
                f"<div class='card'>\n{svg_content}\n</div>\n</body>\n</html>\n"
            )
            out_html.write_text(html_content, encoding="utf-8")
        plotted += 1

    return plotted


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--clusters", required=True, type=Path)
    parser.add_argument("--gbks", nargs="*", type=Path, default=[])
    parser.add_argument("--cargo-profiles", nargs="*", type=Path, default=[])
    parser.add_argument("--max-elements", type=int, default=15)
    parser.add_argument("--outdir", required=True, type=Path)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    n = plot_clusters(
        clusters_tsv=args.clusters,
        gbk_slices=args.gbks,
        cargo_profiles=args.cargo_profiles,
        out_dir=args.outdir,
        max_elements=args.max_elements,
    )
    print(f"Generated synteny plots for {n} cluster(s)", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
