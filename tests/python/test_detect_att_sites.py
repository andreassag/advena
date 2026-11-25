"""Unit tests for bin/detect_att_sites.py."""

from __future__ import annotations

import csv
from pathlib import Path

import detect_att_sites

REPO_ROOT = Path(__file__).resolve().parents[2]


def test_detect_att_sites_on_mini_genome_cargo(tmp_path: Path) -> None:
    """Verify 16 bp direct repeat attL/attR and trnS-Lys tRNA hotspot are detected."""
    cargo_gb = REPO_ROOT / "tests" / "fixtures" / "mini_genome_cargo.gbff"
    me_tsv = tmp_path / "mini_cargo_detected_ME.tsv"
    me_tsv.write_text(
        "ICE_IME_id\tGenome_accession\tIntegrase_upstream\tIntegrase_downstream\t"
        "Start_of_most_upstream_SP\tStop_of_most_downstream_SP\n"
        "ID_1\tNC_000964_CARGO\tBSU_04800\t-\t505\t14764\n",
        encoding="utf-8",
    )

    out_tsv = tmp_path / "att_sites.tsv"
    rc = detect_att_sites.main(
        [
            "--sample",
            "mini_cargo",
            "--detected-me",
            str(me_tsv),
            "--genbank",
            str(cargo_gb),
            "--flank-bp",
            "1000",
            "--no-blast",
            "-o",
            str(out_tsv),
        ]
    )
    assert rc == 0

    with open(out_tsv, encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    assert len(rows) == 1
    r = rows[0]
    assert r["me_id"] == "ID_1"
    assert r["integrase_side"] == "upstream"
    assert int(r["att_repeat_length"]) >= 16
    assert "ATTCCAGTCCACGAGC" in r["att_repeat_seq"]
    assert "tRNA:trnS-Lys" in r["att_target_feature"]
    assert int(r["att_target_distance_bp"]) == 0


def test_detect_att_sites_skips_elements_without_integrase(tmp_path: Path) -> None:
    cargo_gb = REPO_ROOT / "tests" / "fixtures" / "mini_genome_cargo.gbff"
    me_tsv = tmp_path / "no_int_ME.tsv"
    me_tsv.write_text(
        "ICE_IME_id\tGenome_accession\tIntegrase_upstream\tIntegrase_downstream\t"
        "Start_of_most_upstream_SP\tStop_of_most_downstream_SP\n"
        "ID_3\tNC_000964_CARGO\t-\t-\t505\t14764\n",
        encoding="utf-8",
    )
    rows = detect_att_sites.detect_att_sites(
        sample="s1",
        detected_me_tsv=me_tsv,
        genbank_path=cargo_gb,
        use_blastn=False,
    )
    assert rows == []


def test_direct_repeat_orientation_ignores_inverted_only() -> None:
    """Verify direct repeats are detected and inverted-only pairs are ignored."""
    direct = "ACGTACGTCCGGTTAA"
    # Reverse complement of direct is TTAACCGGACGTACGT (not equal to direct)
    revcomp = "TTAACCGGACGTACGT"
    left = "A" * 40 + direct + "C" * 40
    right_inv = "G" * 40 + revcomp + "T" * 40
    right_dir = "G" * 40 + direct + "T" * 40

    hits_inv = detect_att_sites.find_direct_repeats_kmer(left, 1, right_inv, 500, min_len=12)
    assert not any(direct in h.seq for h in hits_inv)

    hits_dir = detect_att_sites.find_direct_repeats_kmer(left, 1, right_dir, 500, min_len=12)
    assert any(h.seq == direct and h.length == 16 for h in hits_dir)


def test_direct_repeat_rejects_homopolymer_and_merges_gbk_targets(tmp_path: Path) -> None:
    """Reject low-complexity homopolymer tracts and merge GFF + GenBank targets."""
    poly_t = "TTTCTTTTATTTTTTTATTTTT"
    hits = detect_att_sites.find_direct_repeats_kmer(
        "ACGT" * 5 + poly_t + "CGTA" * 5,
        1,
        "GGCC" * 5 + poly_t + "TTAA" * 5,
        500,
        min_len=12,
    )
    assert not any(h.seq == poly_t for h in hits)

    # Provide a GFF lacking tRNA features; tRNA:trnS-Lys from GenBank must still be found
    cargo_gb = REPO_ROOT / "tests" / "fixtures" / "mini_genome_cargo.gbff"
    me_tsv = tmp_path / "mini_cargo_detected_ME.tsv"
    me_tsv.write_text(
        "ICE_IME_id\tGenome_accession\tIntegrase_upstream\tIntegrase_downstream\t"
        "Start_of_most_upstream_SP\tStop_of_most_downstream_SP\n"
        "ID_1\tNC_000964_CARGO\tBSU_04800\t-\t505\t14764\n",
        encoding="utf-8",
    )
    gff_no_trna = tmp_path / "icescreen_only.gff"
    gff_no_trna.write_text(
        "##gff-version 3\nNC_000964_CARGO\tICEscreen\tCDS\t505\t1611\t.\t-\t0\tgene=rpsI\n",
        encoding="utf-8",
    )
    rows = detect_att_sites.detect_att_sites(
        sample="mini_cargo",
        detected_me_tsv=me_tsv,
        genbank_path=cargo_gb,
        gff_path=gff_no_trna,
        flank_bp=1000,
        use_blastn=False,
    )
    assert len(rows) == 1
    assert "tRNA:trnS-Lys" in rows[0]["att_target_feature"]
    assert int(rows[0]["att_target_distance_bp"]) == 0
