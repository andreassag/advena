"""Unit tests for bin/extract_me_sequences.py."""

from __future__ import annotations

import csv
import gzip
from pathlib import Path

import extract_me_sequences
from gbk_utils import parse_genbank, read_fasta

REPO_ROOT = Path(__file__).resolve().parents[2]


def _write_mini_me_tsv(path: Path, contig: str = "NC_000964.3", start: int = 505, end: int = 12229):
    cols = [
        "ICE_IME_id",
        "Segment_number",
        "Genome_accession",
        "Category_of_element",
        "Start_of_most_upstream_SP",
        "Stop_of_most_downstream_SP",
    ]
    with open(path, "w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t")
        w.writeheader()
        w.writerow(
            {
                "ICE_IME_id": "ID_1",
                "Segment_number": "1",
                "Genome_accession": contig,
                "Category_of_element": "Complete ICE (I+C+R+V)",
                "Start_of_most_upstream_SP": str(start),
                "Stop_of_most_downstream_SP": str(end),
            }
        )


def _write_mini_sp_tsv(path: Path) -> None:
    path.write_text(
        "ICE_IME_id\tCDS_locus_tag\tCDS_Protein_type\n"
        "ID_1\tBSU_04800\tTyrosine integrase\n"
        "ID_1\tBSU_04860\tCoupling protein\n"
        "ID_1\tBSU_04870\tRelaxase\n"
        "ID_1\tBSU_04940\tVirB4\n",
        encoding="utf-8",
    )


def test_extract_me_sequences_default_and_flank_clamping(tmp_path: Path) -> None:
    mini_gb = REPO_ROOT / "tests" / "fixtures" / "mini_genome.gbff"
    me_tsv = tmp_path / "mini_detected_ME.tsv"
    sp_tsv = tmp_path / "mini_detected_SP_withMEIds.tsv"
    _write_mini_me_tsv(me_tsv, start=505, end=12229)
    _write_mini_sp_tsv(sp_tsv)

    # Default flank = 0
    out0 = tmp_path / "out_flank0"
    n = extract_me_sequences.extract_me_sequences(
        sample="mini",
        prefix="mini",
        detected_me_tsv=me_tsv,
        sp_with_me_tsv=sp_tsv,
        genbank_path=mini_gb,
        out_dir=out0,
        flank_bp=0,
    )
    assert n == 1
    assert (out0 / "mini_me_count.txt").read_text(encoding="utf-8").strip() == "1"

    bed_line = (out0 / "mini_ME_regions.bed").read_text(encoding="utf-8").strip().split("\t")
    assert bed_line[0] == "NC_000964.3"
    assert int(bed_line[1]) == 504  # 0-based
    assert int(bed_line[2]) == 12229
    assert bed_line[3] == "ID_1"

    fna = read_fasta(out0 / "mini_ME_regions.fna.gz")
    assert "mini__ID_1" in fna
    assert len(fna["mini__ID_1"]) == 12229 - 505 + 1

    faa = read_fasta(out0 / "mini_ME_cds.faa.gz")
    assert "BSU_04800" in faa
    assert "BSU_04940" in faa
    assert "BSU_04950" not in faa  # outside element bounds

    with open(out0 / "mini_ME_cds_manifest.tsv", encoding="utf-8") as fh:
        manifest = list(csv.DictReader(fh, delimiter="\t"))
    by_lt = {r["locus_tag"]: r for r in manifest}
    assert by_lt["BSU_04800"]["is_signature_protein"] == "TRUE"
    assert by_lt["BSU_04800"]["sp_type"] == "Tyrosine integrase"
    assert by_lt["BSU_04810"]["is_signature_protein"] == "FALSE"

    # Sliced GenBank parses cleanly and preserves ICEscreen SP qualifier
    slice_gb = out0 / "slices" / "mini__ID_1.gbk.gz"
    assert slice_gb.exists()
    sliced_recs = parse_genbank(slice_gb)
    assert len(sliced_recs) == 1
    assert len(sliced_recs[0].sequence) == 12229 - 505 + 1
    sp_notes = [
        n_txt
        for f in sliced_recs[0].features
        for n_txt in f.qualifiers.get("note", [])
        if n_txt.startswith("ICEscreen_SP:")
    ]
    assert len(sp_notes) == 4

    # Large flank clamps cleanly to [1, contig_len=13000]
    out_flank = tmp_path / "out_flank5000"
    extract_me_sequences.extract_me_sequences(
        sample="mini",
        prefix="mini",
        detected_me_tsv=me_tsv,
        sp_with_me_tsv=sp_tsv,
        genbank_path=mini_gb,
        out_dir=out_flank,
        flank_bp=5000,
    )
    bed_clamped = (
        (out_flank / "mini_ME_regions.bed").read_text(encoding="utf-8").strip().split("\t")
    )
    assert int(bed_clamped[1]) == 0
    assert int(bed_clamped[2]) == 13000


def test_extract_me_sequences_zero_elements(tmp_path: Path) -> None:
    mini_gb = REPO_ROOT / "tests" / "fixtures" / "mini_genome.gbff"
    empty_me = tmp_path / "zero_detected_ME.tsv"
    empty_me.write_text(
        "ICE_IME_id\tGenome_accession\tStart_of_most_upstream_SP\n", encoding="utf-8"
    )

    out_dir = tmp_path / "out_zero"
    n = extract_me_sequences.extract_me_sequences(
        sample="zero",
        prefix="zero",
        detected_me_tsv=empty_me,
        sp_with_me_tsv=None,
        genbank_path=mini_gb,
        out_dir=out_dir,
    )
    assert n == 0
    assert (out_dir / "zero_me_count.txt").read_text(encoding="utf-8").strip() == "0"
    assert (out_dir / "zero_ME_regions.bed").read_text(encoding="utf-8") == ""
    assert read_fasta(out_dir / "zero_ME_cds.faa.gz") == {}
    with gzip.open(out_dir / "zero_ME_cds.gff3.gz", "rt", encoding="utf-8") as fh:
        assert fh.read().strip() == "##gff-version 3"


def test_extract_me_sequences_genome_scope(tmp_path: Path) -> None:
    mini_gb = REPO_ROOT / "tests" / "fixtures" / "mini_genome.gbff"
    me_tsv = tmp_path / "mini_detected_ME.tsv"
    _write_mini_me_tsv(me_tsv, start=505, end=12229)

    out_dir = tmp_path / "out_genome"
    extract_me_sequences.extract_me_sequences(
        sample="mini",
        prefix="mini",
        detected_me_tsv=me_tsv,
        sp_with_me_tsv=None,
        genbank_path=mini_gb,
        out_dir=out_dir,
        scope="genome",
    )
    with open(out_dir / "mini_ME_cds_manifest.tsv", encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    by_lt = {r["locus_tag"]: r for r in rows}
    assert by_lt["BSU_04800"]["me_id"] == "ID_1"
    assert by_lt["BSU_04950"]["me_id"] == "GENOME_BG"

    # In genome scope, cargo_query.fna.gz must carry full contig IDs matching the GFF3,
    # while ME_regions.fna.gz remains sliced to extracted mobile elements only.
    cargo_fna = read_fasta(out_dir / "mini_cargo_query.fna.gz")
    me_fna = read_fasta(out_dir / "mini_ME_regions.fna.gz")
    assert "NC_000964.3" in cargo_fna
    assert len(cargo_fna["NC_000964.3"]) == 13000
    assert "mini__ID_1" in me_fna
    with gzip.open(out_dir / "mini_ME_cds.gff3.gz", "rt", encoding="utf-8") as fh:
        gff_contigs = {ln.split("\t")[0] for ln in fh if ln.strip() and not ln.startswith("#")}
    assert gff_contigs == {"NC_000964.3"}


def test_gbk_utils_join_complement_and_escaped_quotes(tmp_path: Path) -> None:
    """Verify join(complement(...)) strand/splicing and escaped quotes in qualifiers."""
    gb = tmp_path / "edge.gbff"
    gb.write_text(
        "LOCUS       TEST_EDGE                30 bp    DNA     linear   BCT 01-JAN-2026\n"
        "ACCESSION   TEST_EDGE\n"
        "VERSION     TEST_EDGE.1\n"
        "FEATURES             Location/Qualifiers\n"
        "     CDS             join(complement(17..22),complement(7..12))\n"
        '                     /locus_tag="EDGE_1"\n'
        '                     /product="protein ""quoted"""\n'
        "ORIGIN\n"
        "        1 aaaaaaatta aaaaaattac aaaaaaaaaa\n"
        "//\n",
        encoding="utf-8",
    )
    recs = parse_genbank(gb)
    assert len(recs) == 1
    feat = recs[0].features[0]
    assert feat.strand == "-"
    assert feat.get_first("product") == 'protein "quoted"'
    from gbk_utils import extract_feature_dna

    # 17..22 = TTACAA -> rc = TTGTAA; 7..12 = ATTAAA -> rc = TTTAAT
    dna = extract_feature_dna(recs[0], feat)
    assert len(dna) == 12
    assert dna == "TTGTAATTTAAT"
