"""Unit tests for bin/contextualize_me.py."""

from __future__ import annotations

import csv
from pathlib import Path

import contextualize_me


def test_contextualize_me_flags_plasmid_prophage_is_integron_att(tmp_path: Path) -> None:
    me_tsv = tmp_path / "s1_detected_ME.tsv"
    me_tsv.write_text(
        "ICE_IME_id\tGenome_accession\tCategory_of_element\tStart_of_most_upstream_SP\tStop_of_most_downstream_SP\n"
        "ID_1\tNC_000964.3\tComplete ICE (I+C+R+V)\t10000\t30000\n"
        "ID_2\tNZ_PLASMID1.1\tConjugation module (C+V+R)\t5000\t20000\n",
        encoding="utf-8",
    )

    mob = tmp_path / "contig_report.txt"
    mob.write_text(
        "sample_id\tmolecule_type\tprimary_cluster_id\tcontig_id\n"
        "s1\tchromosome\t-\tNC_000964.3\n"
        "s1\tplasmid\tAA001\tNZ_PLASMID1.1\n",
        encoding="utf-8",
    )
    g_plas = tmp_path / "plasmid_summary.tsv"
    g_plas.write_text(
        "seq_name\tplasmid_score\nNC_000964.3\t0.05\nNZ_PLASMID1.1\t0.92\n",
        encoding="utf-8",
    )
    g_vir = tmp_path / "virus_summary.tsv"
    g_vir.write_text(
        "seq_name\tlength\tcoordinates\nNC_000964.3|provirus_25000_45000\t20001\t25000-45000\n",
        encoding="utf-8",
    )
    isescan = tmp_path / "isescan.tsv"
    isescan.write_text(
        "seqID\tfamily\tisBegin\tisEnd\n"
        "NC_000964.3\tIS256\t7000\t8300\n"
        "NC_000964.3\tIS3\t15000\t16200\n"
        "NC_000964.3\tIS110\t90000\t91500\n",
        encoding="utf-8",
    )
    integron = tmp_path / "s1.integrons"
    integron.write_text(
        "# integron_finder\n"
        "ID_integron\tID_replicon\tpos_beg\tpos_end\n"
        "integron_01\tNC_000964.3\t18000\t21000\n",
        encoding="utf-8",
    )
    att = tmp_path / "s1_att_sites.tsv"
    att.write_text(
        "me_id\tatt_target_feature\tatt_repeat_length\tatt_repeat_seq\n"
        "ID_1\ttRNA:trnS-Lys\t16\tATTCCAGTCCACGAGC\n",
        encoding="utf-8",
    )

    out_tsv = tmp_path / "s1_detected_ME_contextualized.tsv"
    rc = contextualize_me.main(
        [
            "--detected-me",
            str(me_tsv),
            "--mobsuite",
            str(mob),
            "--genomad-plasmid",
            str(g_plas),
            "--genomad-virus",
            str(g_vir),
            "--isescan",
            str(isescan),
            "--integronfinder",
            str(integron),
            "--att-sites",
            str(att),
            "-o",
            str(out_tsv),
        ]
    )
    assert rc == 0

    with open(out_tsv, encoding="utf-8") as fh:
        rows = {r["ICE_IME_id"]: r for r in csv.DictReader(fh, delimiter="\t")}

    id1 = rows["ID_1"]
    assert id1["plasmid_borne_warning"] == "FALSE"
    assert id1["mobsuite_molecule_type"] == "chromosome"
    assert "provirus_25000_45000" in id1["composite_prophage_ids"]
    assert "IS256" in id1["flanking_is_families"]
    assert "IS3" in id1["flanking_is_families"]
    assert "IS110" not in id1["flanking_is_families"]
    assert id1["internal_is_count"] == "1"
    assert id1["integron_ids"] == "integron_01"
    assert id1["att_target_feature"] == "tRNA:trnS-Lys"
    assert id1["att_repeat_length"] == "16"

    id2 = rows["ID_2"]
    assert id2["plasmid_borne_warning"] == "TRUE"
    assert id2["mobsuite_cluster_id"] == "AA001"
    assert id2["composite_prophage_ids"] == "-"
