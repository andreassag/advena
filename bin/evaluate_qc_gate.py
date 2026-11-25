#!/usr/bin/env python3
"""Evaluate CheckM2, QUAST, and GTDB-Tk outputs and emit unified QC gate tables."""

from __future__ import annotations

import argparse
import csv
import re
import sys
from pathlib import Path

BACILLOTA_RE = re.compile(r"p__(Bacillota|Firmicutes)(?:_[A-Z])?(?:;|$)", re.IGNORECASE)


def _strip_ext(name: str) -> str:
    stem = Path(name).name
    for suffix in (".fasta.gz", ".fna.gz", ".fa.gz", ".gbff.gz", ".gbk.gz", ".gb.gz"):
        if stem.endswith(suffix):
            return stem[: -len(suffix)]
    for suffix in (".fasta", ".fna", ".fa", ".gbff", ".gbk", ".gb"):
        if stem.endswith(suffix):
            return stem[: -len(suffix)]
    return stem


def _parse_float(value: str | None) -> float | None:
    if value is None:
        return None
    val = value.strip()
    if not val or val in {"-", "NA", "N/A", "None"}:
        return None
    try:
        return float(val)
    except ValueError:
        return None


def parse_checkm2(paths: list[Path]) -> dict[str, dict[str, float | None]]:
    """Parse one or more CheckM2 quality_report.tsv files."""
    metrics: dict[str, dict[str, float | None]] = {}
    for path in paths:
        if not path.exists() or path.stat().st_size == 0:
            continue
        with open(path, encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            for row in reader:
                raw_name = row.get("Name") or row.get("Bin Id") or ""
                if not raw_name:
                    continue
                sample = _strip_ext(raw_name)
                metrics[sample] = {
                    "completeness": _parse_float(row.get("Completeness")),
                    "contamination": _parse_float(row.get("Contamination")),
                }
    return metrics


def parse_quast(paths: list[Path]) -> dict[str, dict[str, float | None]]:
    """Parse one or more QUAST transposed_report.tsv or report.tsv files."""
    metrics: dict[str, dict[str, float | None]] = {}
    for path in paths:
        if not path.exists() or path.stat().st_size == 0:
            continue
        with open(path, encoding="utf-8", newline="") as handle:
            rows = list(csv.reader(handle, delimiter="\t"))
        if len(rows) < 2:
            continue
        header = rows[0]
        if header and header[0].strip().lower() == "assembly":
            # transposed_report.tsv: rows are assemblies
            col_idx = {col.strip(): i for i, col in enumerate(header)}
            for row in rows[1:]:
                if not row:
                    continue
                sample = _strip_ext(row[0])
                metrics[sample] = {
                    "n_contigs": _parse_float(
                        row[col_idx["# contigs"]] if "# contigs" in col_idx else None
                    ),
                    "total_length": _parse_float(
                        row[col_idx["Total length"]] if "Total length" in col_idx else None
                    ),
                    "n50": _parse_float(row[col_idx["N50"]] if "N50" in col_idx else None),
                }
        else:
            # standard report.tsv: columns are assemblies
            assemblies = [_strip_ext(col) for col in header[1:]]
            by_metric: dict[str, list[str]] = {r[0].strip(): r[1:] for r in rows[1:] if r}
            for idx, sample in enumerate(assemblies):
                metrics[sample] = {
                    "n_contigs": _parse_float(
                        by_metric.get("# contigs", [None] * len(assemblies))[idx]
                    ),
                    "total_length": _parse_float(
                        by_metric.get("Total length", [None] * len(assemblies))[idx]
                    ),
                    "n50": _parse_float(by_metric.get("N50", [None] * len(assemblies))[idx]),
                }
    return metrics


def parse_gtdbtk(paths: list[Path]) -> dict[str, str]:
    """Parse one or more GTDB-Tk summary TSV files."""
    taxonomy: dict[str, str] = {}
    for path in paths:
        if not path.exists() or path.stat().st_size == 0:
            continue
        with open(path, encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            for row in reader:
                genome = row.get("user_genome") or ""
                classification = row.get("classification") or ""
                if genome:
                    taxonomy[_strip_ext(genome)] = classification.strip()
    return taxonomy


def is_bacillota(classification: str) -> bool | None:
    """Return True if GTDB-Tk classification is Bacillota/Firmicutes, False if other, None if NA."""
    if not classification or classification in {"-", "Unclassified", "N/A"}:
        return None
    return bool(BACILLOTA_RE.search(classification))


def evaluate_qc(
    samples: list[str],
    checkm2_paths: list[Path],
    quast_paths: list[Path],
    gtdbtk_paths: list[Path],
    min_completeness: float = 90.0,
    max_contamination: float = 5.0,
    strict_qc: bool = False,
    filter_non_bacillota: bool = False,
) -> tuple[list[dict[str, str]], list[str]]:
    """Evaluate all samples and return (summary_rows, passed_sample_ids)."""
    cm2 = parse_checkm2(checkm2_paths)
    quast = parse_quast(quast_paths)
    gtdb = parse_gtdbtk(gtdbtk_paths)

    all_samples = list(dict.fromkeys(samples))
    for extra in list(cm2) + list(quast) + list(gtdb):
        if extra not in all_samples:
            all_samples.append(extra)

    summary_rows: list[dict[str, str]] = []
    passed: list[str] = []

    for sample in all_samples:
        c_info = cm2.get(sample, {})
        q_info = quast.get(sample, {})
        tax = gtdb.get(sample, "")

        comp = c_info.get("completeness")
        cont = c_info.get("contamination")
        n_contigs = q_info.get("n_contigs")
        total_len = q_info.get("total_length")
        n50 = q_info.get("n50")
        bacillota_flag = is_bacillota(tax)

        flags: list[str] = []
        excluded = False

        if comp is not None and comp < min_completeness:
            prefix = "FAIL" if strict_qc else "WARN"
            flags.append(f"{prefix}_LOW_COMPLETENESS")
            if strict_qc:
                excluded = True

        if cont is not None and cont > max_contamination:
            prefix = "FAIL" if strict_qc else "WARN"
            flags.append(f"{prefix}_HIGH_CONTAMINATION")
            if strict_qc:
                excluded = True

        if bacillota_flag is False:
            prefix = "FAIL" if filter_non_bacillota else "WARN"
            flags.append(f"{prefix}_NON_BACILLOTA")
            if filter_non_bacillota:
                excluded = True

        status = ";".join(flags) if flags else "PASS"
        if not excluded:
            passed.append(sample)

        summary_rows.append(
            {
                "sample": sample,
                "completeness": f"{comp:.2f}" if comp is not None else "NA",
                "contamination": f"{cont:.2f}" if cont is not None else "NA",
                "n_contigs": str(int(n_contigs)) if n_contigs is not None else "NA",
                "total_length": str(int(total_len)) if total_len is not None else "NA",
                "n50": str(int(n50)) if n50 is not None else "NA",
                "gtdb_classification": tax or "NA",
                "is_bacillota": (
                    "TRUE"
                    if bacillota_flag is True
                    else ("FALSE" if bacillota_flag is False else "NA")
                ),
                "qc_status": status,
                "included_downstream": "FALSE" if excluded else "TRUE",
            }
        )

    return summary_rows, passed


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--samples", nargs="*", default=[], help="Expected sample IDs")
    parser.add_argument("--checkm2", nargs="*", type=Path, default=[], help="CheckM2 report TSVs")
    parser.add_argument("--quast", nargs="*", type=Path, default=[], help="QUAST report TSVs")
    parser.add_argument("--gtdbtk", nargs="*", type=Path, default=[], help="GTDB-Tk summary TSVs")
    parser.add_argument("--min-completeness", type=float, default=90.0)
    parser.add_argument("--max-contamination", type=float, default=5.0)
    parser.add_argument("--strict-qc", action="store_true")
    parser.add_argument("--filter-non-bacillota", action="store_true")
    parser.add_argument("--out-summary", required=True, type=Path)
    parser.add_argument("--out-passed", required=True, type=Path)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    rows, passed = evaluate_qc(
        samples=args.samples,
        checkm2_paths=args.checkm2,
        quast_paths=args.quast,
        gtdbtk_paths=args.gtdbtk,
        min_completeness=args.min_completeness,
        max_contamination=args.max_contamination,
        strict_qc=args.strict_qc,
        filter_non_bacillota=args.filter_non_bacillota,
    )

    fieldnames = [
        "sample",
        "completeness",
        "contamination",
        "n_contigs",
        "total_length",
        "n50",
        "gtdb_classification",
        "is_bacillota",
        "qc_status",
        "included_downstream",
    ]
    args.out_summary.parent.mkdir(parents=True, exist_ok=True)
    with open(args.out_summary, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    args.out_passed.parent.mkdir(parents=True, exist_ok=True)
    args.out_passed.write_text("\n".join(passed) + ("\n" if passed else ""), encoding="utf-8")

    for row in rows:
        if row["qc_status"] != "PASS":
            print(
                f"[QC GATE] sample={row['sample']} status={row['qc_status']} "
                f"included={row['included_downstream']}",
                file=sys.stderr,
            )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
