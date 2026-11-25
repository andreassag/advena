"""Unit tests for bin/intersect_bgc_me.py and bin/merge_cargo_annotations.py."""

from __future__ import annotations

import csv
import gzip
import json
from pathlib import Path

import intersect_bgc_me
import merge_cargo_annotations


def test_intersect_bgc_and_merge_cargo_annotations(tmp_path: Path) -> None:
    me_tsv = tmp_path / "s1_detected_ME.tsv"
    me_tsv.write_text(
        "ICE_IME_id\tGenome_accession\tStart_of_most_upstream_SP\tStop_of_most_downstream_SP\n"
        "ID_1\tNC_000964.3\t500\t12000\n",
        encoding="utf-8",
    )

    manifest = tmp_path / "s1_ME_cds_manifest.tsv"
    manifest.write_text(
        "\t".join(
            [
                "sample",
                "me_id",
                "me_category",
                "contig",
                "locus_tag",
                "protein_id",
                "orig_start",
                "orig_end",
                "strand",
                "is_signature_protein",
                "sp_type",
                "product",
            ]
        )
        + "\n"
        + "s1\tID_1\tComplete ICE\tNC_000964.3\tLOC_1\tWP_1\t500\t1600\t-\t"
        + "TRUE\tTyrosine integrase\tintegrase\n"
        + "s1\tID_1\tComplete ICE\tNC_000964.3\tLOC_2\tWP_2\t5000\t6500\t+\tFALSE\t-\ttetM\n"
        + "s1\tID_1\tComplete ICE\tNC_000964.3\tLOC_3\tWP_3\t8000\t9500\t+\tFALSE\t-\tbgc\n",
        encoding="utf-8",
    )

    # antiSMASH JSON with one overlapping bacteriocin BGC and one distant non-overlapping BGC
    anti_json = tmp_path / "s1.json"
    anti_json.write_text(
        json.dumps(
            {
                "records": [
                    {
                        "id": "NC_000964.3",
                        "areas": [
                            {"start": 7500, "end": 10000, "products": ["lanthipeptide-class-i"]},
                            {"start": 50000, "end": 60000, "products": ["NRPS"]},
                        ],
                    }
                ]
            }
        ),
        encoding="utf-8",
    )

    bgc_tsv = tmp_path / "s1_bgc_me_overlaps.tsv"
    assert (
        intersect_bgc_me.main(
            [
                "--sample",
                "s1",
                "--detected-me",
                str(me_tsv),
                "--manifest",
                str(manifest),
                "--antismash-json",
                str(anti_json),
                "-o",
                str(bgc_tsv),
            ]
        )
        == 0
    )

    with open(bgc_tsv, encoding="utf-8") as fh:
        bgc_rows = list(csv.DictReader(fh, delimiter="\t"))
    assert len(bgc_rows) == 1
    assert bgc_rows[0]["bgc_products"] == "lanthipeptide-class-i"
    assert bgc_rows[0]["overlapping_locus_tags"] == "LOC_3"

    # Build AMR, eggNOG, DefenseFinder, and ICEscreen GFF3 inputs
    amr_tsv = tmp_path / "amr.tsv"
    amr_tsv.write_text(
        "sample\tme_id\tcontig\tlocus_tag\tcanonical_symbol\telement_type\t"
        "consensus_drug_class\tcard_aro_accession\tdetected_by\n"
        "s1\tID_1\tNC_000964.3\tLOC_2\ttet(M)\tAMR\ttetracycline\tARO:3000186\tboth\n",
        encoding="utf-8",
    )
    eggnog_tsv = tmp_path / "eggnog.tsv"
    eggnog_tsv.write_text(
        "sample\tme_id\tcontig\tlocus_tag\tbacteria_og\tbacillota_og\t"
        "narrowest_og\tcog_category\tdescription\tcargo_role\n"
        "s1\tID_1\tNC_000964.3\tLOC_2\tCOG0480@2\t2IQ7R@1239\t"
        "2IQ7R@1239\tJ\tTetM\tCARGO_CHARACTERIZED_OG\n",
        encoding="utf-8",
    )
    def_sys = tmp_path / "defense_systems.tsv"
    def_sys.write_text("sys_id\ttype\tsubtype\nSYS_1\tRM\tRM_Type_II\n", encoding="utf-8")
    def_genes = tmp_path / "defense_genes.tsv"
    def_genes.write_text(
        "hit_id\tsys_id\tgene_name\nLOC_3\tSYS_1\tType_II_MTases\n", encoding="utf-8"
    )

    ice_gff = tmp_path / "s1_icescreen.gff.gz"
    with gzip.open(ice_gff, "wt", encoding="utf-8") as fh:
        fh.write(
            "##gff-version 3\n"
            "NC_000964.3\tRefSeq\tCDS\t5000\t6500\t.\t+\t0\tID=cds2;locus_tag=LOC_2\n"
            "NC_000964.3\tRefSeq\tCDS\t8000\t9500\t.\t+\t0\tID=cds3;locus_tag=LOC_3\n"
        )

    out_prof = tmp_path / "s1_me_cargo_profile.tsv"
    out_gff = tmp_path / "s1_advena.gff.gz"
    assert (
        merge_cargo_annotations.main(
            [
                "--manifest",
                str(manifest),
                "--icescreen-gff",
                str(ice_gff),
                "--amr",
                str(amr_tsv),
                "--eggnog",
                str(eggnog_tsv),
                "--defense-genes",
                str(def_genes),
                "--defense-systems",
                str(def_sys),
                "--bgc",
                str(bgc_tsv),
                "--out-profile",
                str(out_prof),
                "--out-gff",
                str(out_gff),
            ]
        )
        == 0
    )

    with open(out_prof, encoding="utf-8") as fh:
        prof = {r["locus_tag"]: r for r in csv.DictReader(fh, delimiter="\t")}
    assert prof["LOC_2"]["amr_gene"] == "tet(M)"
    assert prof["LOC_2"]["card_aro"] == "ARO:3000186"
    assert prof["LOC_2"]["eggnog_bacillota_og"] == "2IQ7R@1239"
    assert prof["LOC_3"]["defense_system"] == "RM_Type_II"
    assert prof["LOC_3"]["bgc_type"] == "lanthipeptide-class-i"

    with gzip.open(out_gff, "rt", encoding="utf-8") as fh:
        gff_text = fh.read()
    assert "amr_gene=tet(M)" in gff_text
    assert "card_aro=ARO:3000186" in gff_text
    assert "bacillota_og=2IQ7R@1239" in gff_text
    assert "defense_system=RM_Type_II" in gff_text
    assert "bgc_type=lanthipeptide-class-i" in gff_text
