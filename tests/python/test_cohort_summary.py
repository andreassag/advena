"""Unit tests for bin/cluster_mes.py, bin/build_cohort_summary.py, and bin/plot_me_synteny.py."""

from __future__ import annotations

import csv
import gzip
from pathlib import Path

import build_cohort_summary
import cluster_mes
import extract_me_sequences
import plot_me_synteny

REPO_ROOT = Path(__file__).resolve().parents[2]


def test_cluster_mes_skani_and_short_ime_mmseqs_fallback(tmp_path: Path) -> None:
    bed1 = tmp_path / "s1_ME_regions.bed"
    bed1.write_text(
        "NC_1\t0\t20000\tID_1\t0\t+\tComplete ICE\nNC_1\t30000\t33000\tID_2\t0\t+\tComplete IME\n",
        encoding="utf-8",
    )
    bed2 = tmp_path / "s2_ME_regions.bed"
    bed2.write_text(
        "NC_2\t0\t20500\tID_1\t0\t+\tComplete ICE\nNC_2\t30000\t33200\tID_2\t0\t+\tComplete IME\n",
        encoding="utf-8",
    )

    m1 = tmp_path / "s1_ME_cds_manifest.tsv"
    m1.write_text(
        "sample\tme_id\tme_category\tlocus_tag\n"
        "s1\tID_2\tComplete IME\tS1_P1\n"
        "s1\tID_2\tComplete IME\tS1_P2\n",
        encoding="utf-8",
    )
    m2 = tmp_path / "s2_ME_cds_manifest.tsv"
    m2.write_text(
        "sample\tme_id\tme_category\tlocus_tag\n"
        "s2\tID_2\tComplete IME\tS2_P1\n"
        "s2\tID_2\tComplete IME\tS2_P2\n",
        encoding="utf-8",
    )

    # skani clusters the 20 kb ICEs (ANI 98.5%, AF >= 90%)
    skani = tmp_path / "skani.tsv"
    skani.write_text(
        "Ref_file\tQuery_file\tANI\tAlign_fraction_ref\tAlign_fraction_query\n"
        "s1__ID_1.fna\ts2__ID_1.fna\t98.5\t94.0\t91.2\n",
        encoding="utf-8",
    )
    # MMseqs2 protein clusters link the 3 kb short IMEs (Jaccard = 1.0)
    mmseqs = tmp_path / "mmseqs_cluster.tsv"
    mmseqs.write_text(
        "S1_P1\tS1_P1\nS1_P1\tS2_P1\nS1_P2\tS1_P2\nS1_P2\tS2_P2\n",
        encoding="utf-8",
    )

    out_tsv = tmp_path / "me_clusters.tsv"
    out_json = tmp_path / "cluster_members.json"
    rc = cluster_mes.main(
        [
            "--beds",
            str(bed1),
            str(bed2),
            "--manifests",
            str(m1),
            str(m2),
            "--skani",
            str(skani),
            "--mmseqs",
            str(mmseqs),
            "--out-tsv",
            str(out_tsv),
            "--out-json",
            str(out_json),
        ]
    )
    assert rc == 0

    with open(out_tsv, encoding="utf-8") as fh:
        rows = {r["element_key"]: r for r in csv.DictReader(fh, delimiter="\t")}
    assert rows["s1__ID_1"]["me_cluster_id"] == rows["s2__ID_1"]["me_cluster_id"]
    assert rows["s1__ID_1"]["clustering_evidence"] == "skani_ani"
    assert rows["s1__ID_2"]["me_cluster_id"] == rows["s2__ID_2"]["me_cluster_id"]
    assert rows["s1__ID_2"]["clustering_evidence"] == "mmseqs_jaccard"
    assert rows["s1__ID_1"]["me_cluster_id"] != rows["s1__ID_2"]["me_cluster_id"]


def test_cluster_mes_handles_identical_locus_tags_with_sample_prefix(tmp_path: Path) -> None:
    """Two samples sharing identical locus_tag IDs cluster cleanly via sample|locus_tag prefix."""
    bed1 = tmp_path / "mini_cargo_ME_regions.bed"
    bed1.write_text("NC_1\t0\t3500\tID_1\t0\t+\tComplete IME\n", encoding="utf-8")
    bed2 = tmp_path / "mini_ref_ME_regions.bed"
    bed2.write_text("NC_1\t0\t3400\tID_1\t0\t+\tComplete IME\n", encoding="utf-8")

    m1 = tmp_path / "mini_cargo_ME_cds_manifest.tsv"
    m1.write_text(
        "sample\tme_id\tme_category\tlocus_tag\n"
        "mini_cargo\tID_1\tComplete IME\tBSU_04800\n"
        "mini_cargo\tID_1\tComplete IME\tBSU_04810\n",
        encoding="utf-8",
    )
    m2 = tmp_path / "mini_ref_ME_cds_manifest.tsv"
    m2.write_text(
        "sample\tme_id\tme_category\tlocus_tag\n"
        "mini_ref\tID_1\tComplete IME\tBSU_04800\n"
        "mini_ref\tID_1\tComplete IME\tBSU_04810\n",
        encoding="utf-8",
    )

    mmseqs = tmp_path / "pan_me_cluster.tsv"
    mmseqs.write_text(
        "mini_cargo|BSU_04800\tmini_cargo|BSU_04800\n"
        "mini_cargo|BSU_04800\tmini_ref|BSU_04800\n"
        "mini_cargo|BSU_04810\tmini_cargo|BSU_04810\n"
        "mini_cargo|BSU_04810\tmini_ref|BSU_04810\n",
        encoding="utf-8",
    )

    rows, members = cluster_mes.cluster_elements(
        bed_paths=[bed1, bed2],
        manifest_paths=[m1, m2],
        mmseqs_path=mmseqs,
    )
    assert len(members) == 1
    assert sorted(next(iter(members.values()))) == ["mini_cargo__ID_1", "mini_ref__ID_1"]
    assert all(r["clustering_evidence"] == "mmseqs_jaccard" for r in rows)


def test_build_cohort_summary_and_pan_og_matrix(tmp_path: Path) -> None:
    me_hdr = (
        "ICE_IME_id\tGenome_accession\tCategory_of_element\tICEline_format\t"
        "ICE_consensus_superfamily_SP_conj_module\tICE_consensus_family_SP_conj_module\t"
        "Start_of_most_upstream_SP\tStop_of_most_downstream_SP\tplasmid_borne_warning\t"
        "composite_prophage_ids\tflanking_is_families\tatt_target_feature\n"
    )
    me1 = tmp_path / "s1_detected_ME_contextualized.tsv"
    me1.write_text(
        me_hdr + "ID_1\tNC_1\tComplete ICE (I+C+R+V)\ti7.C0.R6.V\tTn916\tICEBs1\t"
        "505\t12229\tFALSE\t-\tIS256\ttRNA:trnS-Lys\n",
        encoding="utf-8",
    )
    me2 = tmp_path / "s2_detected_ME_contextualized.tsv"
    me2.write_text(
        me_hdr + "ID_1\tNC_2\tComplete ICE (I+C+R+V)\ti7.C0.R6.V\tTn916\tICEBs1\t"
        "505\t12229\tTRUE\t-\t-\t-\n",
        encoding="utf-8",
    )

    cargo_hdr = (
        "sample\tme_id\tlocus_tag\tamr_gene\tdefense_system\tbgc_type\t"
        "eggnog_bacteria_og\teggnog_bacillota_og\tcog_category\teggnog_description\tcargo_role\n"
    )
    c1 = tmp_path / "s1_me_cargo_profile.tsv"
    c1.write_text(
        cargo_hdr + "s1\tID_1\tL1\t-\t-\t-\tCOG0582@2\t2IGK4@1239\tL\tIntegrase\tICE_BACKBONE_SP\n"
        "s1\tID_1\tL2\ttet(M)\tRM_Type_II\t-\tCOG0480@2\t2IQ7R@1239\t"
        "J\tTetM\tCARGO_CHARACTERIZED_OG\n",
        encoding="utf-8",
    )
    c2 = tmp_path / "s2_me_cargo_profile.tsv"
    c2.write_text(
        cargo_hdr + "s2\tID_1\tM1\t-\t-\t-\tCOG0582@2\t2IGK4@1239\tL\tIntegrase\tICE_BACKBONE_SP\n",
        encoding="utf-8",
    )

    out_dir = tmp_path / "cohort"
    rc = build_cohort_summary.main(
        [
            "--contextualized-me",
            str(me1),
            str(me2),
            "--cargo-profiles",
            str(c1),
            str(c2),
            "--outdir",
            str(out_dir),
        ]
    )
    assert rc == 0

    with open(out_dir / "cohort_me_cargo_summary.tsv", encoding="utf-8") as fh:
        summary = {r["sample"]: r for r in csv.DictReader(fh, delimiter="\t")}
    assert len(summary) == 2
    assert summary["s1"]["amr_genes"] == "tet(M)"
    assert summary["s1"]["defense_systems"] == "RM_Type_II"
    assert summary["s2"]["plasmid_borne_warning"] == "TRUE"

    with open(out_dir / "og_presence_absence.tsv", encoding="utf-8") as fh:
        ogs = {r["og_id"]: r for r in csv.DictReader(fh, delimiter="\t")}
    assert ogs["2IGK4@1239"]["classification"] == "CORE_BACKBONE"
    assert ogs["2IGK4@1239"]["s1|ID_1"] == "1"
    assert ogs["2IGK4@1239"]["s2|ID_1"] == "1"
    assert ogs["2IQ7R@1239"]["classification"] == "SHELL_CARGO"
    assert ogs["2IQ7R@1239"]["s1|ID_1"] == "1"
    assert ogs["2IQ7R@1239"]["s2|ID_1"] == "0"

    assert (out_dir / "icescreen_summary_mqc.yaml").exists()
    assert "tet(M)" in (out_dir / "amr_card_heatmap_mqc.tsv").read_text(encoding="utf-8")
    assert "advena_cog_distribution" in (out_dir / "cog_distribution_mqc.tsv").read_text(
        encoding="utf-8"
    )
    assert "Plasmid_borne_MEs" in (out_dir / "mge_context_mqc.tsv").read_text(encoding="utf-8")


def test_plot_me_synteny_renders_svg_and_html(tmp_path: Path) -> None:
    # Slice mini_genome and mini_genome_cargo to create two real GenBank slices
    mini_gb = REPO_ROOT / "tests" / "fixtures" / "mini_genome.gbff"
    cargo_gb = REPO_ROOT / "tests" / "fixtures" / "mini_genome_cargo.gbff"

    me_cols = (
        "ICE_IME_id\tGenome_accession\tCategory_of_element\t"
        "Start_of_most_upstream_SP\tStop_of_most_downstream_SP\n"
    )
    me1 = tmp_path / "m1.tsv"
    me1.write_text(me_cols + "ID_1\tNC_000964.3\tComplete ICE\t505\t12229\n", encoding="utf-8")
    me2 = tmp_path / "m2.tsv"
    me2.write_text(me_cols + "ID_1\tNC_000964_CARGO\tComplete ICE\t505\t14764\n", encoding="utf-8")
    extract_me_sequences.extract_me_sequences("ref", "ref", me1, None, mini_gb, tmp_path / "e1")
    extract_me_sequences.extract_me_sequences(
        "cargo", "cargo", me2, None, cargo_gb, tmp_path / "e2"
    )

    s1_gbk = tmp_path / "e1" / "slices" / "ref__ID_1.gbk.gz"
    s2_gbk = tmp_path / "e2" / "slices" / "cargo__ID_1.gbk.gz"
    assert s1_gbk.exists() and s2_gbk.exists()
    with gzip.open(s2_gbk, "rt", encoding="utf-8") as fh:
        assert "CARGO_TETM" in fh.read()

    clusters = tmp_path / "me_clusters.tsv"
    clusters.write_text(
        "me_cluster_id\telement_key\tsample\tme_id\n"
        "CLUSTER_001\tref__ID_1\tref\tID_1\n"
        "CLUSTER_001\tcargo__ID_1\tcargo\tID_1\n",
        encoding="utf-8",
    )
    cargo_prof = tmp_path / "cargo_profile.tsv"
    cargo_prof.write_text(
        "sample\tlocus_tag\tis_signature_protein\tsp_type\tamr_gene\t"
        "defense_system\tbgc_type\tcargo_role\teggnog_bacillota_og\n"
        "cargo\tBSU_04800\tTRUE\tIntTyr\t-\t-\t-\tICE_BACKBONE_SP\t-\n"
        "cargo\tCARGO_TETM\tFALSE\t-\ttet(M)\t-\t-\tCARGO_CHARACTERIZED_OG\t2IQ7R@1239\n",
        encoding="utf-8",
    )

    out_dir = tmp_path / "synteny"
    rc = plot_me_synteny.main(
        [
            "--clusters",
            str(clusters),
            "--gbks",
            str(s1_gbk),
            str(s2_gbk),
            "--cargo-profiles",
            str(cargo_prof),
            "--outdir",
            str(out_dir),
        ]
    )
    assert rc == 0

    svg_path = out_dir / "CLUSTER_001_synteny.svg"
    html_path = out_dir / "CLUSTER_001_synteny.html"
    assert svg_path.exists() and html_path.exists()
    svg_text = svg_path.read_text(encoding="utf-8")
    assert "<svg " in svg_text
    assert "CLUSTER_001" in svg_text
    assert "tet(M)" in svg_text
    assert "<!DOCTYPE html>" in html_path.read_text(encoding="utf-8")
