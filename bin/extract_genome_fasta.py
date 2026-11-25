#!/usr/bin/env python3
"""Extract nucleotide FASTA records from a GenBank (.gb/.gbk/.gbff[.gz]) file."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from gbk_utils import open_text, parse_genbank, wrap_fasta


def extract_fasta(genbank_path: Path, output_path: Path) -> int:
    """Write nucleotide FASTA records from GenBank; return number of records written."""
    records = parse_genbank(genbank_path)
    if not records:
        raise ValueError(f"No GenBank records found in {genbank_path}")

    written = 0
    with open_text(output_path, "wt") as out:
        for rec in records:
            if not rec.sequence:
                continue
            header = f">{rec.record_id}"
            if rec.definition:
                header += f" {rec.definition}"
            out.write(header + "\n")
            out.write(wrap_fasta(rec.sequence) + "\n")
            written += 1

    if written == 0:
        raise ValueError(f"No nucleotide sequences found in {genbank_path}")
    return written


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-i", "--genbank", required=True, type=Path, help="Input GenBank file")
    parser.add_argument("-o", "--output", required=True, type=Path, help="Output FASTA (.fna/.gz)")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    count = extract_fasta(args.genbank, args.output)
    print(f"Extracted {count} sequence(s) to {args.output}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
