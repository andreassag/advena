#!/usr/bin/env python3
"""Shared pure-Python GenBank, FASTA, and GFF3 utilities for advena helper scripts."""

from __future__ import annotations

import gzip
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import IO, TextIO

_COMPLEMENT_TABLE = str.maketrans(
    "ACGTacgtNnRrYySsWwKkMmBbDdHhVv", "TGCAtgcaNnYyRrSsWwMmKkVvHhDdBb"
)

# NCBI Translation Table 11 (Bacterial, Archaeal and Plant Plastid Code)
_CODON_TABLE_11: dict[str, str] = {
    "TTT": "F",
    "TTC": "F",
    "TTA": "L",
    "TTG": "L",
    "TCT": "S",
    "TCC": "S",
    "TCA": "S",
    "TCG": "S",
    "TAT": "Y",
    "TAC": "Y",
    "TAA": "*",
    "TAG": "*",
    "TGT": "C",
    "TGC": "C",
    "TGA": "*",
    "TGG": "W",
    "CTT": "L",
    "CTC": "L",
    "CTA": "L",
    "CTG": "L",
    "CCT": "P",
    "CCC": "P",
    "CCA": "P",
    "CCG": "P",
    "CAT": "H",
    "CAC": "H",
    "CAA": "Q",
    "CAG": "Q",
    "CGT": "R",
    "CGC": "R",
    "CGA": "R",
    "CGG": "R",
    "ATT": "I",
    "ATC": "I",
    "ATA": "I",
    "ATG": "M",
    "ACT": "T",
    "ACC": "T",
    "ACA": "T",
    "ACG": "T",
    "AAT": "N",
    "AAC": "N",
    "AAA": "K",
    "AAG": "K",
    "AGT": "S",
    "AGC": "S",
    "AGA": "R",
    "AGG": "R",
    "GTT": "V",
    "GTC": "V",
    "GTA": "V",
    "GTG": "V",
    "GCT": "A",
    "GCC": "A",
    "GCA": "A",
    "GCG": "A",
    "GAT": "D",
    "GAC": "D",
    "GAA": "E",
    "GAG": "E",
    "GGT": "G",
    "GGC": "G",
    "GGA": "G",
    "GGG": "G",
}


def open_text(path: str | Path, mode: str = "rt") -> TextIO | IO[str]:
    """Open a plain or gzipped text file transparently."""
    path_str = str(path)
    if path_str.endswith(".gz"):
        return gzip.open(path_str, mode, encoding="utf-8")
    return open(path_str, mode, encoding="utf-8")


def reverse_complement(seq: str) -> str:
    """Return the reverse complement of a DNA sequence."""
    return seq.translate(_COMPLEMENT_TABLE)[::-1]


def translate_bacterial(dna_seq: str) -> str:
    """Translate a bacterial CDS (NCBI table 11), trimming trailing stop codon."""
    seq = dna_seq.upper().replace("U", "T")
    aa: list[str] = []
    for i in range(0, len(seq) - 2, 3):
        codon = seq[i : i + 3]
        residue = _CODON_TABLE_11.get(codon, "X")
        if i == 0 and codon in {"TTG", "CTG", "ATT", "ATC", "ATA", "ATG", "GTG"}:
            residue = "M"
        aa.append(residue)
    protein = "".join(aa)
    if protein.endswith("*"):
        protein = protein[:-1]
    return protein


@dataclass
class GBFeature:
    """A single feature entry from a GenBank record."""

    feature_type: str
    location: str
    start: int  # 1-based inclusive
    end: int  # 1-based inclusive
    strand: str  # '+' or '-'
    parts: list[tuple[int, int]] = field(default_factory=list)
    qualifiers: dict[str, list[str]] = field(default_factory=dict)

    def get_first(self, key: str, default: str = "") -> str:
        vals = self.qualifiers.get(key)
        return vals[0] if vals else default


@dataclass
class GBRecord:
    """A single GenBank sequence record."""

    record_id: str
    locus: str
    accession: str
    version: str
    definition: str
    organism: str
    topology: str
    sequence: str
    features: list[GBFeature] = field(default_factory=list)


def parse_location(loc_str: str) -> tuple[int, int, str, list[tuple[int, int]]]:
    """Parse a GenBank feature location string into (start, end, strand, parts)."""
    cleaned = re.sub(r"\s+", "", loc_str)
    strand = "-" if "complement(" in cleaned else "+"
    spans: list[tuple[int, int]] = []
    for m in re.finditer(r"[<>]?(\d+)(?:\.\.[<>]?(\d+))?", cleaned):
        s = int(m.group(1))
        e = int(m.group(2)) if m.group(2) else s
        if s > e:
            s, e = e, s
        spans.append((s, e))
    if not spans:
        return (1, 1, strand, [(1, 1)])
    return (min(p[0] for p in spans), max(p[1] for p in spans), strand, spans)


def _finalize_qualifier(
    qualifiers: dict[str, list[str]], key: str | None, value_lines: list[str]
) -> None:
    if key is None:
        return
    if not value_lines:
        qualifiers.setdefault(key, []).append("true")
        return
    raw = "\n".join(value_lines)
    if raw.startswith('"') and raw.endswith('"') and len(raw) >= 2:
        inner = raw[1:-1].replace('""', '"')
        inner = "".join(inner.split()) if key == "translation" else re.sub(r"\s*\n\s*", " ", inner)
        qualifiers.setdefault(key, []).append(inner)
    else:
        qualifiers.setdefault(key, []).append(raw.strip())


def parse_genbank(path: str | Path) -> list[GBRecord]:
    """Parse all GenBank records from a plain or .gz GenBank file."""
    records: list[GBRecord] = []
    with open_text(path, "rt") as handle:
        locus = ""
        accession = ""
        version = ""
        definition_lines: list[str] = []
        organism = ""
        topology = "linear"
        features: list[GBFeature] = []
        seq_chunks: list[str] = []
        section = "HEADER"
        in_definition = False

        cur_feat_type: str | None = None
        cur_loc_lines: list[str] = []
        cur_quals: dict[str, list[str]] = {}
        cur_qual_key: str | None = None
        cur_qual_val_lines: list[str] = []
        in_qual_quote = False

        def flush_feature() -> None:
            nonlocal cur_feat_type, cur_loc_lines, cur_quals
            nonlocal cur_qual_key, cur_qual_val_lines, in_qual_quote
            if cur_feat_type is not None:
                _finalize_qualifier(cur_quals, cur_qual_key, cur_qual_val_lines)
                loc_str = "".join(cur_loc_lines)
                start, end, strand, parts = parse_location(loc_str)
                features.append(
                    GBFeature(
                        feature_type=cur_feat_type,
                        location=loc_str,
                        start=start,
                        end=end,
                        strand=strand,
                        parts=parts,
                        qualifiers=cur_quals,
                    )
                )
            cur_feat_type = None
            cur_loc_lines = []
            cur_quals = {}
            cur_qual_key = None
            cur_qual_val_lines = []
            in_qual_quote = False

        for raw_line in handle:
            line = raw_line.rstrip("\r\n")
            if line == "//":
                flush_feature()
                if locus or version or accession or seq_chunks:
                    rec_id = version or accession or locus or f"contig_{len(records) + 1}"
                    records.append(
                        GBRecord(
                            record_id=rec_id,
                            locus=locus or rec_id,
                            accession=accession or rec_id,
                            version=version or accession or locus or rec_id,
                            definition=" ".join(definition_lines).strip(),
                            organism=organism.strip(),
                            topology=topology,
                            sequence="".join(seq_chunks).upper(),
                            features=features,
                        )
                    )
                locus = ""
                accession = ""
                version = ""
                definition_lines = []
                organism = ""
                topology = "linear"
                features = []
                seq_chunks = []
                section = "HEADER"
                in_definition = False
                continue

            if line.startswith("ORIGIN"):
                flush_feature()
                section = "ORIGIN"
                continue

            if section == "ORIGIN":
                seq_part = re.sub(r"[^A-Za-z]", "", line)
                if seq_part:
                    seq_chunks.append(seq_part)
                continue

            if line.startswith("FEATURES"):
                section = "FEATURES"
                in_definition = False
                continue

            if section == "HEADER":
                if line.startswith("LOCUS"):
                    parts = line.split()
                    if len(parts) >= 2:
                        locus = parts[1]
                    if "circular" in parts:
                        topology = "circular"
                    in_definition = False
                elif line.startswith("DEFINITION"):
                    definition_lines = [line[12:].strip()]
                    in_definition = True
                elif line.startswith("ACCESSION"):
                    parts = line[12:].split()
                    if parts:
                        accession = parts[0]
                    in_definition = False
                elif line.startswith("VERSION"):
                    parts = line[12:].split()
                    if parts:
                        version = parts[0]
                    in_definition = False
                elif line.startswith("  ORGANISM"):
                    organism = line[12:].strip()
                    in_definition = False
                elif in_definition and line.startswith("            "):
                    definition_lines.append(line.strip())
                else:
                    in_definition = False
                continue

            if section == "FEATURES":
                if len(line) > 5 and line.startswith("     ") and not line.startswith("      "):
                    flush_feature()
                    tokens = line.strip().split(None, 1)
                    cur_feat_type = tokens[0]
                    cur_loc_lines = [tokens[1].strip()] if len(tokens) > 1 else []
                    continue

                if cur_feat_type is None:
                    continue

                content = line[21:] if len(line) >= 21 else line.strip()
                if in_qual_quote:
                    cur_qual_val_lines.append(content)
                    trailing_quotes = len(content) - len(content.rstrip('"'))
                    if trailing_quotes % 2 == 1:
                        in_qual_quote = False
                    continue

                if content.startswith("/"):
                    _finalize_qualifier(cur_quals, cur_qual_key, cur_qual_val_lines)
                    cur_qual_val_lines = []
                    if "=" in content:
                        qkey, qval = content[1:].split("=", 1)
                        cur_qual_key = qkey.strip()
                        cur_qual_val_lines = [qval]
                        if qval.startswith('"'):
                            inner_after_open = qval[1:]
                            trailing_quotes = len(inner_after_open) - len(
                                inner_after_open.rstrip('"')
                            )
                            if trailing_quotes % 2 == 0:
                                in_qual_quote = True
                    else:
                        cur_qual_key = content[1:].strip()
                        cur_qual_val_lines = []
                elif cur_qual_key is None:
                    cur_loc_lines.append(content.strip())
                else:
                    cur_qual_val_lines.append(content)

        # Flush trailing record if file lacked final '//'
        flush_feature()
        if locus or version or accession or seq_chunks:
            rec_id = version or accession or locus or f"contig_{len(records) + 1}"
            records.append(
                GBRecord(
                    record_id=rec_id,
                    locus=locus or rec_id,
                    accession=accession or rec_id,
                    version=version or accession or locus or rec_id,
                    definition=" ".join(definition_lines).strip(),
                    organism=organism.strip(),
                    topology=topology,
                    sequence="".join(seq_chunks).upper(),
                    features=features,
                )
            )
    return records


def extract_feature_dna(record: GBRecord, feat: GBFeature) -> str:
    """Extract spliced DNA sequence for a feature (reverse-complemented on '-' strand)."""
    cleaned_loc = re.sub(r"\s+", "", feat.location)
    if feat.strand == "-" and cleaned_loc.startswith("join("):
        # join(complement(b1..b2),complement(a1..a2)) specifies pieces in 5'->3' order
        return "".join(
            reverse_complement(record.sequence[max(0, s - 1) : min(len(record.sequence), e)])
            for s, e in feat.parts
        )
    pieces = [record.sequence[max(0, s - 1) : min(len(record.sequence), e)] for s, e in feat.parts]
    seq = "".join(pieces)
    if feat.strand == "-":
        return reverse_complement(seq)
    return seq


def extract_cds_protein(record: GBRecord, feat: GBFeature) -> str:
    """Return protein sequence from /translation qualifier or bacterial translation fallback."""
    trans = feat.get_first("translation", "").strip()
    if trans:
        return trans
    dna = extract_feature_dna(record, feat)
    if len(dna) >= 3:
        return translate_bacterial(dna)
    return ""


def wrap_fasta(seq: str, width: int = 60) -> str:
    """Wrap sequence string to fixed line width."""
    return "\n".join(seq[i : i + width] for i in range(0, len(seq), width))


def _format_qualifier(key: str, value: str) -> list[str]:
    if value == "true" and key in {"pseudo", "ribosomal_slippage"}:
        return [f"                     /{key}"]
    escaped = value.replace('"', '""')
    text = f'/{key}="{escaped}"'
    max_len = 58
    if len(text) <= max_len:
        return [f"                     {text}"]
    lines: list[str] = []
    pos = 0
    while pos < len(text):
        chunk = text[pos : pos + max_len]
        lines.append(f"                     {chunk}")
        pos += max_len
    return lines


def format_genbank_record(
    locus_name: str,
    definition: str,
    organism: str,
    sequence: str,
    features: list[GBFeature],
    topology: str = "linear",
) -> str:
    """Serialize a GenBank record into standard flat-file format."""
    clean_locus = re.sub(r"\s+", "_", locus_name)[:20]
    seq_lower = sequence.lower()
    locus_line = (
        f"LOCUS       {clean_locus:<20} {len(seq_lower):>7} bp    "
        f"DNA     {topology:<8} BCT 01-JAN-2026"
    )
    out: list[str] = [
        locus_line,
        f"DEFINITION  {definition or clean_locus}.",
        f"ACCESSION   {clean_locus}",
        f"VERSION     {clean_locus}",
        "KEYWORDS    .",
        f"SOURCE      {organism or 'Unknown'}",
        f"  ORGANISM  {organism or 'Unknown'}",
        "            Bacteria.",
        "FEATURES             Location/Qualifiers",
        f"     source          1..{max(1, len(seq_lower))}",
        f'                     /organism="{organism or "Unknown"}"',
    ]

    for feat in features:
        loc = (
            f"complement({feat.start}..{feat.end})"
            if feat.strand == "-"
            else f"{feat.start}..{feat.end}"
        )
        out.append(f"     {feat.feature_type:<16}{loc}")
        for qkey, qvals in feat.qualifiers.items():
            for qval in qvals:
                out.extend(_format_qualifier(qkey, qval))

    out.append("ORIGIN")
    for idx in range(0, len(seq_lower), 60):
        chunk = seq_lower[idx : idx + 60]
        blocks = " ".join(chunk[j : j + 10] for j in range(0, len(chunk), 10))
        out.append(f"{idx + 1:>9} {blocks}")
    out.append("//")
    return "\n".join(out) + "\n"


def read_fasta(path: str | Path) -> dict[str, str]:
    """Read a plain or gzipped FASTA file into an ordered {id: sequence} dict."""
    seqs: dict[str, list[str]] = {}
    cur_id: str | None = None
    with open_text(path, "rt") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                cur_id = line[1:].split()[0]
                seqs[cur_id] = []
            elif cur_id is not None:
                seqs[cur_id].append(line)
    return {k: "".join(v) for k, v in seqs.items()}


def parse_gff3_attributes(attr_str: str) -> dict[str, str]:
    """Parse a GFF3 column-9 attributes string into a dict."""
    attrs: dict[str, str] = {}
    for part in attr_str.strip().split(";"):
        if not part or "=" not in part:
            continue
        k, v = part.split("=", 1)
        attrs[k.strip()] = v.strip()
    return attrs


def format_gff3_attributes(attrs: dict[str, str]) -> str:
    """Serialize a dict into a GFF3 column-9 attributes string."""
    return ";".join(f"{k}={v}" for k, v in attrs.items() if v != "")
