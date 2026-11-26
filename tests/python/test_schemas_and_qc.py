"""Unit tests for JSON schema validation, FASTA extraction, and QC gate evaluation."""

from __future__ import annotations

import csv
import gzip
import subprocess
from pathlib import Path

import pytest

import evaluate_qc_gate
import extract_genome_fasta

REPO_ROOT = Path(__file__).resolve().parents[2]


def test_schemas_valid_and_synced() -> None:
    res = subprocess.run(
        ["python3", str(REPO_ROOT / "tests" / "python" / "validate_schemas.py")],
        capture_output=True,
        text=True,
        check=False,
    )
    assert res.returncode == 0, f"validate_schemas.py failed:\n{res.stdout}\n{res.stderr}"


def test_editorconfig_compliance() -> None:
    """Verify tracked files match .editorconfig rules (respecting .editorconfig-checker.json)."""
    tracked = subprocess.check_output(["git", "ls-files"], cwd=REPO_ROOT, text=True).splitlines()
    excludes = ("assets/", "modules/nf-core/", "tests/fixtures/")
    errors: list[str] = []

    for rel in tracked:
        if any(rel.startswith(ex) for ex in excludes):
            continue
        p = REPO_ROOT / rel
        if not p.is_file():
            continue
        raw = p.read_bytes()
        if b"\r\n" in raw:
            errors.append(f"{rel}: CRLF line endings")
        try:
            text = raw.decode("utf-8")
        except UnicodeDecodeError:
            errors.append(f"{rel}: non-UTF-8 file not excluded")
            continue
        if text and not text.endswith("\n"):
            errors.append(f"{rel}: missing final newline")
        lines = text.splitlines()
        ext = p.suffix
        if p.name != "Makefile" and ext not in {".md", ".diff"}:
            for i, line in enumerate(lines, 1):
                if line != line.rstrip(" \t"):
                    errors.append(f"{rel}:{i}: trailing whitespace")
        if p.name == "Makefile":
            for i, line in enumerate(lines, 1):
                if line.startswith(" "):
                    errors.append(f"{rel}:{i}: space indent in Makefile")
        else:
            for i, line in enumerate(lines, 1):
                ws = line[: len(line) - len(line.lstrip())]
                if "\t" in ws:
                    errors.append(f"{rel}:{i}: tab indent")
        if ext == ".py":
            for i, line in enumerate(lines, 1):
                if len(line) > 120:
                    errors.append(f"{rel}:{i}: line length {len(line)} > 120")
                stripped = line.lstrip(" ")
                indent = len(line) - len(stripped)
                if stripped and indent % 4 != 0:
                    errors.append(f"{rel}:{i}: Python indent {indent} not multiple of 4")
        elif ext in {".json", ".yml", ".yaml", ".toml", ".md"}:
            for i, line in enumerate(lines, 1):
                stripped = line.lstrip(" ")
                indent = len(line) - len(stripped)
                if stripped and indent % 2 != 0:
                    errors.append(f"{rel}:{i}: {ext} indent {indent} not multiple of 2")

    assert not errors, "EditorConfig violations:\n" + "\n".join(errors[:25])


def test_extract_genome_fasta_plain_and_gz(tmp_path: Path) -> None:
    mini_gb = REPO_ROOT / "tests" / "fixtures" / "mini_genome.gbff"
    out_fna = tmp_path / "mini.fna"
    n = extract_genome_fasta.extract_fasta(mini_gb, out_fna)
    assert n == 1
    lines = out_fna.read_text(encoding="utf-8").splitlines()
    assert lines[0].startswith(">NC_000964.3")
    seq = "".join(lines[1:])
    assert len(seq) == 13000

    # Gzipped input & output roundtrip
    gz_in = tmp_path / "mini.gbff.gz"
    with gzip.open(gz_in, "wt", encoding="utf-8") as fh:
        fh.write(mini_gb.read_text(encoding="utf-8"))
    out_gz = tmp_path / "mini.fna.gz"
    assert extract_genome_fasta.extract_fasta(gz_in, out_gz) == 1
    with gzip.open(out_gz, "rt", encoding="utf-8") as fh:
        assert fh.readline().startswith(">NC_000964.3")


def test_extract_genome_fasta_empty_raises(tmp_path: Path) -> None:
    empty_gb = tmp_path / "empty.gbff"
    empty_gb.write_text("", encoding="utf-8")
    with pytest.raises(ValueError, match="No GenBank records"):
        extract_genome_fasta.extract_fasta(empty_gb, tmp_path / "out.fna")


def test_evaluate_qc_gate_warn_vs_strict_and_phylum(tmp_path: Path) -> None:
    cm2 = tmp_path / "quality_report.tsv"
    cm2.write_text(
        "Name\tCompleteness\tContamination\n"
        "s_pass.fna\t99.5\t0.8\n"
        "s_lowcomp.fna\t82.0\t1.0\n"
        "s_contam.fna\t98.0\t7.5\n"
        "s_ecoli.fna\t99.9\t0.1\n",
        encoding="utf-8",
    )
    quast = tmp_path / "transposed_report.tsv"
    quast.write_text(
        "Assembly\t# contigs\tTotal length\tN50\n"
        "s_pass\t1\t4200000\t4200000\n"
        "s_lowcomp\t250\t2100000\t15000\n",
        encoding="utf-8",
    )
    gtdb = tmp_path / "gtdbtk.bac120.summary.tsv"
    gtdb.write_text(
        "user_genome\tclassification\n"
        "s_pass\td__Bacteria;p__Bacillota;c__Bacilli\n"
        "s_lowcomp\td__Bacteria;p__Firmicutes_A;c__Clostridia\n"
        "s_ecoli\td__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria\n",
        encoding="utf-8",
    )

    samples = ["s_pass", "s_lowcomp", "s_contam", "s_ecoli"]

    # Default mode: warn only, all samples included downstream
    rows_warn, passed_warn = evaluate_qc_gate.evaluate_qc(
        samples=samples,
        checkm2_paths=[cm2],
        quast_paths=[quast],
        gtdbtk_paths=[gtdb],
        strict_qc=False,
        filter_non_bacillota=False,
    )
    assert passed_warn == samples
    by_s = {r["sample"]: r for r in rows_warn}
    assert by_s["s_pass"]["qc_status"] == "PASS"
    assert by_s["s_pass"]["is_bacillota"] == "TRUE"
    assert by_s["s_lowcomp"]["qc_status"] == "WARN_LOW_COMPLETENESS"
    assert by_s["s_lowcomp"]["is_bacillota"] == "TRUE"  # Firmicutes_A recognized
    assert by_s["s_contam"]["qc_status"] == "WARN_HIGH_CONTAMINATION"
    assert by_s["s_ecoli"]["qc_status"] == "WARN_NON_BACILLOTA"
    assert by_s["s_ecoli"]["is_bacillota"] == "FALSE"

    # Strict QC + phylum filter mode: only s_pass survives
    rows_strict, passed_strict = evaluate_qc_gate.evaluate_qc(
        samples=samples,
        checkm2_paths=[cm2],
        quast_paths=[quast],
        gtdbtk_paths=[gtdb],
        strict_qc=True,
        filter_non_bacillota=True,
    )
    assert passed_strict == ["s_pass"]
    strict_map = {r["sample"]: r for r in rows_strict}
    assert strict_map["s_lowcomp"]["included_downstream"] == "FALSE"
    assert strict_map["s_contam"]["included_downstream"] == "FALSE"
    assert strict_map["s_ecoli"]["included_downstream"] == "FALSE"

    # CLI execution writes TSV and passed list cleanly
    out_tsv = tmp_path / "qc_gate_summary.tsv"
    out_passed = tmp_path / "passed.txt"
    rc = evaluate_qc_gate.main(
        [
            "--samples",
            *samples,
            "--checkm2",
            str(cm2),
            "--quast",
            str(quast),
            "--gtdbtk",
            str(gtdb),
            "--strict-qc",
            "--out-summary",
            str(out_tsv),
            "--out-passed",
            str(out_passed),
        ]
    )
    assert rc == 0
    with open(out_tsv, encoding="utf-8") as fh:
        written = list(csv.DictReader(fh, delimiter="\t"))
    assert len(written) == 4
    assert out_passed.read_text(encoding="utf-8").splitlines() == ["s_pass", "s_ecoli"]
