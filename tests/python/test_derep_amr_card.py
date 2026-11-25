"""Unit tests for bin/derep_amr_card.py (AMRFinderPlus --plus + CARD RGI dereplication)."""

from __future__ import annotations

import csv
from pathlib import Path

import derep_amr_card


def _write_manifest(path: Path) -> None:
    path.write_text(
        "sample\tme_id\tme_category\tcontig\tlocus_tag\tprotein_id\t"
        "orig_start\torig_end\tstrand\tis_signature_protein\tsp_type\tproduct\n"
        "s1\tID_1\tComplete ICE\tNC_1\tLOC_TETM\tWP_1\t1000\t2920\t+\tFALSE\t-\ttetM\n"
        "s1\tID_1\tComplete ICE\tNC_1\tLOC_CADA\tWP_2\t3000\t4000\t+\tFALSE\t-\tcadA\n"
        "s1\tID_1\tComplete ICE\tNC_1\tLOC_VANA\tWP_3\t5000\t6030\t-\tFALSE\t-\tvanA\n"
        "s1\tID_1\tComplete ICE\tNC_1\tLOC_ADJ\tWP_4\t6300\t7300\t-\tFALSE\t-\tadj\n",
        encoding="utf-8",
    )


def test_derep_locus_tag_merge_and_stress_preservation(tmp_path: Path) -> None:
    manifest = tmp_path / "manifest.tsv"
    _write_manifest(manifest)

    amr_tsv = tmp_path / "amrfinder.tsv"
    amr_tsv.write_text(
        "Protein id\tContig id\tStart\tStop\tStrand\tElement symbol\tElement name\tScope\t"
        "Type\tSubtype\tClass\tSubclass\tMethod\t% Coverage of reference\t"
        "% Identity to reference\tClosest reference accession\n"
        "LOC_TETM\tNC_1\t1000\t2920\t+\ttet(M)\ttetracycline resistance protein\tcore\t"
        "AMR\tAMR\tTETRACYCLINE\tTETRACYCLINE\tEXACTP\t100.0\t100.0\tWP_0001\n"
        "LOC_CADA\tNC_1\t3000\t4000\t+\tcadA\tcadmium ATPase\tplus\t"
        "STRESS\tMETAL\tCADMIUM\tCADMIUM\tBLASTP\t99.0\t98.5\tWP_0002\n",
        encoding="utf-8",
    )

    rgi_txt = tmp_path / "card_rgi.txt"
    rgi_txt.write_text(
        "ORF_ID\tContig\tStart\tStop\tOrientation\tCut_Off\tBest_Hit_ARO\tBest_Identities\t"
        "ARO\tDrug Class\tResistance Mechanism\tPercentage Length of Reference Sequence\n"
        "LOC_TETM\tNC_1\t1000\t2920\t+\tPerfect\ttet(M)\t100.0\t3000186\t"
        "tetracycline antibiotic\tantibiotic target protection\t100.0\n"
        "LOC_VANA\tNC_1\t5000\t6030\t-\tStrict\tvanA\t97.5\t3002639\t"
        "glycopeptide antibiotic\tantibiotic target alteration\t99.5\n"
        "LOC_ADJ\tNC_1\t6300\t7300\t-\tLoose\tadeF\t45.0\t3000777\t"
        "fluoroquinolone antibiotic\tantibiotic efflux\t60.0\n",
        encoding="utf-8",
    )

    out_tsv = tmp_path / "derep.tsv"
    rc = derep_amr_card.main(
        [
            "--sample",
            "s1",
            "--manifest",
            str(manifest),
            "--amrfinder",
            str(amr_tsv),
            "--card-rgi",
            str(rgi_txt),
            "-o",
            str(out_tsv),
        ]
    )
    assert rc == 0

    with open(out_tsv, encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))

    # Loose hit excluded by default; 3 dereplicated calls remain: tet(M), cadA (STRESS/METAL), vanA
    assert len(rows) == 3
    by_lt = {r["locus_tag"]: r for r in rows}

    tetm = by_lt["LOC_TETM"]
    assert tetm["canonical_symbol"] == "tet(M)"
    assert tetm["detected_by"] == "both"
    assert tetm["card_aro_accession"] == "ARO:3000186"
    assert tetm["resistance_mechanism"] == "antibiotic target protection"
    assert tetm["consensus_drug_class"] == "tetracycline"
    assert tetm["representative_tier"] == "3"

    cada = by_lt["LOC_CADA"]
    assert cada["canonical_symbol"] == "cadA"
    assert cada["element_type"] == "STRESS"
    assert cada["element_subtype"] == "METAL"
    assert cada["detected_by"] == "amrfinderplus"

    vana = by_lt["LOC_VANA"]
    assert vana["canonical_symbol"] == "vanA"
    assert vana["card_aro_accession"] == "ARO:3002639"
    assert vana["detected_by"] == "card_rgi"


def test_derep_coordinate_overlap_threshold_and_strand() -> None:
    # Different locus tags (DNA vs protein) with >=80% reciprocal overlap -> MERGE
    h_amr = derep_amr_card.HitRecord(
        source="amrfinderplus",
        sample="s1",
        me_id="ID_1",
        contig="NC_1",
        locus_tag="L1",
        start=1001,
        end=2000,
        strand="+",
        symbol="erm(B)",
        name="ErmB",
        element_type="AMR",
        element_subtype="AMR",
        drug_classes=["macrolide"],
        method_or_cutoff="HMM",
        pct_identity=85.0,
        pct_coverage=90.0,
        tier=1,
    )
    h_card_ov = derep_amr_card.HitRecord(
        source="card_rgi",
        sample="s1",
        me_id="ID_1",
        contig="NC_1",
        locus_tag="ORF_99",
        start=1101,
        end=2050,
        strand="+",
        symbol="ErmB",
        name="ErmB",
        element_type="AMR",
        element_subtype="AMR",
        drug_classes=["lincosamide", "macrolide"],
        method_or_cutoff="Strict",
        pct_identity=96.0,
        pct_coverage=98.0,
        tier=2,
        aro_accession="ARO:3000375",
        aro_name="ErmB",
        resistance_mechanism="antibiotic target alteration",
    )
    # Opposite-strand overlapping gene -> MUST NOT MERGE
    h_opp_strand = derep_amr_card.HitRecord(
        source="card_rgi",
        sample="s1",
        me_id="ID_1",
        contig="NC_1",
        locus_tag="ORF_OPP",
        start=1001,
        end=2000,
        strand="-",
        symbol="tet(L)",
        name="tet(L)",
        element_type="AMR",
        element_subtype="AMR",
        drug_classes=["tetracycline"],
        method_or_cutoff="Perfect",
        pct_identity=100.0,
        pct_coverage=100.0,
        tier=3,
        aro_accession="ARO:3000185",
        aro_name="tet(L)",
    )
    # Low overlap (<80%, e.g. 50%) on same strand -> MUST NOT MERGE
    h_low_ov = derep_amr_card.HitRecord(
        source="amrfinderplus",
        sample="s1",
        me_id="ID_1",
        contig="NC_1",
        locus_tag="L_LOW",
        start=1501,
        end=2500,
        strand="+",
        symbol="qacA",
        name="QacA",
        element_type="STRESS",
        element_subtype="BIOCIDE",
        drug_classes=["quaternary_ammonium"],
        method_or_cutoff="BLASTP",
        pct_identity=95.0,
        pct_coverage=95.0,
        tier=2,
    )

    merged = derep_amr_card.dereplicate_hits(
        [h_amr, h_card_ov, h_opp_strand, h_low_ov], min_overlap=0.80
    )
    assert len(merged) == 3

    erm_row = next(r for r in merged if r["canonical_symbol"] == "erm(B)")
    assert erm_row["detected_by"] == "both"
    # Tier 2 Strict CARD beats Tier 1 HMM AMRFinder as representative
    assert erm_row["representative_source"] == "card_rgi"
    assert erm_row["representative_tier"] == "2"
    assert erm_row["card_aro_accession"] == "ARO:3000375"
    assert erm_row["consensus_drug_class"] == "lincosamide;macrolide"

    assert any(r["canonical_symbol"] == "tet(L)" and r["strand"] == "-" for r in merged)
    assert any(
        r["canonical_symbol"] == "qacA" and r["element_subtype"] == "BIOCIDE" for r in merged
    )


def test_derep_empty_inputs(tmp_path: Path) -> None:
    out_tsv = tmp_path / "empty_derep.tsv"
    rc = derep_amr_card.main(["--sample", "s0", "-o", str(out_tsv)])
    assert rc == 0
    lines = out_tsv.read_text(encoding="utf-8").splitlines()
    assert len(lines) == 1
    assert "canonical_symbol" in lines[0]


def test_derep_lifts_blastx_slice_coords_via_bed(tmp_path: Path) -> None:
    """DNA-only AMRFinderPlus BLASTX hit on element slice dereplicates with CARD RGI protein hit."""
    manifest = tmp_path / "manifest.tsv"
    _write_manifest(manifest)

    # Element ID_1 starts at 1-based position 501 (0-based 500) on NC_1
    me_bed = tmp_path / "s1_ME_regions.bed"
    me_bed.write_text("NC_1\t500\t8000\tID_1\t0\t+\tComplete ICE\n", encoding="utf-8")

    # LOC_TETM is at orig 1000..2920 -> slice-relative 500..2420 on contig s1__ID_1
    amr_tsv = tmp_path / "amrfinder_blastx.tsv"
    amr_tsv.write_text(
        "Protein id\tContig id\tStart\tStop\tStrand\tElement symbol\tElement name\tScope\t"
        "Type\tSubtype\tClass\tSubclass\tMethod\t% Coverage of reference\t"
        "% Identity to reference\tClosest reference accession\n"
        "NA\ts1__ID_1\t500\t2420\t+\ttet(M)\ttetracycline resistance\tcore\t"
        "AMR\tAMR\tTETRACYCLINE\tTETRACYCLINE\tBLASTX\t99.0\t98.0\tWP_0001\n",
        encoding="utf-8",
    )
    rgi_txt = tmp_path / "card_rgi.txt"
    rgi_txt.write_text(
        "ORF_ID\tContig\tStart\tStop\tOrientation\tCut_Off\tBest_Hit_ARO\tBest_Identities\t"
        "ARO\tDrug Class\tResistance Mechanism\tPercentage Length of Reference Sequence\n"
        "LOC_TETM\t\t\t\t\tPerfect\ttet(M)\t100.0\t3000186\t"
        "tetracycline antibiotic\tantibiotic target protection\t100.0\n",
        encoding="utf-8",
    )

    out_tsv = tmp_path / "derep_lifted.tsv"
    rc = derep_amr_card.main(
        [
            "--sample",
            "s1",
            "--manifest",
            str(manifest),
            "--bed",
            str(me_bed),
            "--amrfinder",
            str(amr_tsv),
            "--card-rgi",
            str(rgi_txt),
            "-o",
            str(out_tsv),
        ]
    )
    assert rc == 0
    with open(out_tsv, encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    assert len(rows) == 1
    assert rows[0]["detected_by"] == "both"
    assert rows[0]["contig"] == "NC_1"
    assert rows[0]["me_id"] == "ID_1"
    assert rows[0]["start"] == "1000"
    assert rows[0]["end"] == "2920"
