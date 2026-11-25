"""Unit tests for bin/parse_eggnog_ogs.py."""

from __future__ import annotations

import csv
from pathlib import Path

import parse_eggnog_ogs


def test_parse_og_hierarchy() -> None:
    parsed = parse_eggnog_ogs.parse_og_hierarchy(
        "COG0582@1|root,COG0582@2|Bacteria,2IGK4@1239|Firmicutes,4S999@91061|Bacilli"
    )
    assert parsed["root_cog"] == "COG0582"
    assert parsed["bacteria_og"] == "COG0582@2"
    assert parsed["bacillota_og"] == "2IGK4@1239"
    assert parsed["narrowest_og"] == "4S999@91061"

    empty = parse_eggnog_ogs.parse_og_hierarchy("-")
    assert all(v == "-" for v in empty.values())


def test_parse_emapper_annotations_roles_and_missing(tmp_path: Path) -> None:
    manifest = tmp_path / "manifest.tsv"
    manifest.write_text(
        "sample\tme_id\tme_category\tcontig\tlocus_tag\tprotein_id\t"
        "orig_start\torig_end\tstrand\tis_signature_protein\tsp_type\tproduct\n"
        "s1\tID_1\tComplete ICE\tNC_1\tBSU_04800\tNP_1\t505\t1611\t-\t"
        "TRUE\tTyrosine integrase\tintegrase\n"
        "s1\tID_1\tComplete ICE\tNC_1\tCARGO_TETM\tWP_1\t7000\t8000\t+\tFALSE\t-\ttetM\n"
        "s1\tID_1\tComplete ICE\tNC_1\tCARGO_HYP\tWP_2\t8100\t8400\t+\tFALSE\t-\thypothetical\n"
        "s1\tID_1\tComplete ICE\tNC_1\tCARGO_NOHIT\tWP_3\t8500\t8800\t+\tFALSE\t-\tunknown\n",
        encoding="utf-8",
    )

    emapper = tmp_path / "s1.emapper.annotations"
    emapper.write_text(
        "## emapper-2.1.12\n"
        "#query\tseed_ortholog\tevalue\tscore\teggNOG_OGs\tmax_annot_lvl\tCOG_category\t"
        "Description\tPreferred_name\tGOs\tEC\tKEGG_ko\tKEGG_Pathway\tKEGG_Module\t"
        "KEGG_Reaction\tKEGG_rclass\tBRITE\tKEGG_TC\tCAZy\tBiGG_Reaction\tPFAMs\n"
        "BSU_04800\t1423.BSU04800\t1e-90\t350\t"
        "COG0582@1|root,COG0582@2|Bacteria,2IGK4@1239|Firmicutes\t1239|Firmicutes\tL\t"
        "Tyrosine integrase\tydcL\tGO:0003677\t-\tko:K04763\t-\t-\t-\t-\t-\t-\t-\t-\t"
        "Phage_integrase\n"
        "CARGO_TETM\t1423.TETM\t1e-120\t420\t"
        "COG0480@1|root,COG0480@2|Bacteria,2IQ7R@1239|Firmicutes\t1239|Firmicutes\tJ\t"
        "Tetracycline resistance TetM\ttetM\t-\t-\tko:K18220\tko01501\t-\t-\t-\t-\t-\t-\t-\t"
        "GTP_EFTU\n"
        "CARGO_HYP\t1423.HYP\t1e-10\t75\t"
        "COG9999@1|root,COG9999@2|Bacteria\t2|Bacteria\tS\tFunction unknown\t"
        "-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\tDUF1234\n",
        encoding="utf-8",
    )

    out_tsv = tmp_path / "cargo_ogs.tsv"
    rc = parse_eggnog_ogs.main(
        [
            "--sample",
            "s1",
            "--emapper",
            str(emapper),
            "--manifest",
            str(manifest),
            "-o",
            str(out_tsv),
        ]
    )
    assert rc == 0

    with open(out_tsv, encoding="utf-8") as fh:
        rows = {r["locus_tag"]: r for r in csv.DictReader(fh, delimiter="\t")}

    assert rows["BSU_04800"]["cargo_role"] == "ICE_BACKBONE_SP"
    assert rows["BSU_04800"]["bacillota_og"] == "2IGK4@1239"

    assert rows["CARGO_TETM"]["cargo_role"] == "CARGO_CHARACTERIZED_OG"
    assert rows["CARGO_TETM"]["bacteria_og"] == "COG0480@2"
    assert rows["CARGO_TETM"]["bacillota_og"] == "2IQ7R@1239"
    assert rows["CARGO_TETM"]["kegg_ko"] == "ko:K18220"

    assert rows["CARGO_HYP"]["cargo_role"] == "CARGO_HYPOTHETICAL_OG"
    assert rows["CARGO_HYP"]["cog_category"] == "S"

    assert rows["CARGO_NOHIT"]["cargo_role"] == "CARGO_HYPOTHETICAL_OG"
    assert rows["CARGO_NOHIT"]["bacteria_og"] == "-"
