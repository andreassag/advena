#!/usr/bin/env python3
"""Build deterministic miniature reference DBs and mini_genome_cargo.gbff for offline tests."""

from __future__ import annotations

import json
import sqlite3
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT / "bin"))

from gbk_utils import (  # noqa: E402
    GBFeature,
    format_genbank_record,
    parse_genbank,
    wrap_fasta,
)

# Representative short synthetic/canonical protein fragments for deterministic miniature fixtures
TETM_PROT = (
    "MKIINIGVLAHVDAGKTTLTESLLYNSGAITELGSVDKGTTRTDNTLLERQRGITIQTGITSFQWENTKVNIIDTPGHMD"
    "FLAEVYRSLSVLDGAILLISAKDGVQAQTRILFHALRKMGIPTIFFINKIDQNGIDLSTVYQDIKEKLSAEIVIKQKVEL"
    "YPNMCVTNFTESEQWDTVIEGNDDLLEKYMSGKSLEALELEQEESIRFHNCSLFPVYHGSAKNNIGIDNLIEVITNKFYS"
)
CADA_PROT = (
    "MKKLVILVTGASGGIGKAIALELAKRGAKVIGTARSDEKAKALAEELGAKVVALDVTDEASVKKAVAEVEEKFGRVDVLV"
    "NNAGITRDALLVRMKEEDWNDVINTNLKGAFNVTKAVLPYMKERKSGRIINISSVVGMTGNIGQANYAAAKAGVIGLTKT"
    "AARELASRGITVNAVAPGFIATDMTDDLAEEKILKQIPLGRFGQPEEVANAVLFLASDDASYITGQTINVNGGMYM"
)
RM_MT_PROT = (
    "MKLKDIKGLPTEQIAKVLDSVGKTYEDFVQKFKKLGVTYEELIKKAGLSVAALKKIYETGKISVNEIAKLAGVSKSTVSR"
    "YINKLKKNGYEVQIDPKTGRIKVLDLFAGSGGLGFGLAKAGFKVIGIEKNPEAAKIAEENAKKLGVDNVEFIQGDAKNLP"
)
RM_RE_PROT = (
    "MSNEQIKQLLAEFKELGIEVKDKNGKIVVPDVKLKDGKIVEVSLKGKTYEELIKELGVKVSVAEAKKVLKEVGKTYEDFV"
    "QEFKKLGVTYEELIKKAGLSVAALKKIYETGKISVNEIAKLAGVSKSTVSRYINKLKKNGYEVQIDPKTGRIKVLDLFAG"
)
ERMB_PROT = (
    "MNKNIKYSQNFLTSEKVLNQIIKQLNLKETDTVYEIGTGKGHLTTKLAKISKQVTSIELDSHLFNLSSEKLKLNTRVTLI"
    "HQDILQFQFPNKQRYKIVGNIPYHLSTQIIKKVVFESRASDIYLIVEEGFYKRTLDIHRTLGLLLHTQVSIQQLLKLPAE"
)
VANA_PROT = (
    "MNRIKVAILFGGCSEEHDVSVKSAIEIAANINKEKYEPLYIGITKSGVWKMCEKPCAEWENDNCYSAVLSPDKKMHGLLV"
    "KKNHEYEINHVDVAFSALHGKSGEDGSIQGLFELSGIPFVGCDIQSSAICMDKSLTYIVAKNAGIATPAFWVINKDDRPV"
)
QACA_PROT = (
    "MQTPLKKWLLIVLVTLGFFVMSLDNTVLNVSLPTIKQELGLSTSQVSWIVDIYSLVLAGFIIPLSAFAEKWGRKKALLTG"
    "FALFTFGSLIAFFAPNFATLLVARFLQGIGAGALMPTTLSMIRVIFENPKERATALAVWSIASSIGAVFGPIIGGALLEQ"
)

ATT_REPEAT = "ATTCCAGTCCACGAGC"  # 16 bp direct repeat (attL / attR)

_AA_TO_CODON = {
    "A": "GCT",
    "C": "TGC",
    "D": "GAT",
    "E": "GAA",
    "F": "TTC",
    "G": "GGT",
    "H": "CAT",
    "I": "ATT",
    "K": "AAA",
    "L": "CTG",
    "M": "ATG",
    "N": "AAC",
    "P": "CCG",
    "Q": "CAG",
    "R": "CGT",
    "S": "TCT",
    "T": "ACT",
    "V": "GTT",
    "W": "TGG",
    "Y": "TAT",
}


def _back_translate(protein: str) -> str:
    return "".join(_AA_TO_CODON.get(aa, "GCT") for aa in protein) + "TAA"


def _build_valid_hmm(
    template_hmm: Path,
    name: str,
    desc: str,
    acc: str | None = None,
) -> str:
    """Clone a syntactically valid HMMER3/f profile with updated NAME/ACC/DESC headers."""
    raw_lines = template_hmm.read_text(encoding="utf-8").splitlines()
    out_lines: list[str] = []
    inserted_acc = False
    for line in raw_lines:
        if line.startswith("NAME  "):
            out_lines.append(f"NAME  {name}")
            if acc:
                out_lines.append(f"ACC   {acc}")
                inserted_acc = True
        elif line.startswith("ACC   "):
            if acc and not inserted_acc:
                out_lines.append(f"ACC   {acc}")
                inserted_acc = True
        elif line.startswith("DESC  "):
            out_lines.append(f"DESC  {desc}")
        else:
            out_lines.append(line)
    return "\n".join(out_lines) + "\n"


def build_mini_genome_cargo(fixtures_dir: Path) -> Path:
    """Create mini_genome_cargo.gbff with ICEBs1 + tet(M) + cadA + RM cassette + attL/attR."""
    source_gb = fixtures_dir / "mini_genome.gbff"
    rec = parse_genbank(source_gb)[0]

    # Insert cargo cassette cleanly at position 6950 (between nicK 5773..6831 and ydcS 7096..7365)
    insert_pos = 6950
    spacer = "ATGCGTACGTAGCTAGCTAG"  # 20 bp
    tetm_dna = _back_translate(TETM_PROT)
    cada_dna = _back_translate(CADA_PROT)
    rm_m_dna = _back_translate(RM_MT_PROT)
    rm_r_dna = _back_translate(RM_RE_PROT)

    cargo_seq = (
        spacer + tetm_dna + spacer + cada_dna + spacer + rm_m_dna + spacer + rm_r_dna + spacer
    )
    delta = len(cargo_seq)

    new_seq_list = list(rec.sequence[: insert_pos - 1] + cargo_seq + rec.sequence[insert_pos - 1 :])

    # Place identical 16 bp attL at 3' end of trnS-Lys (147..162) and attR downstream of ICE
    attl_pos = 147
    attr_pos = 12650 + delta
    new_seq_list[attl_pos - 1 : attl_pos - 1 + len(ATT_REPEAT)] = list(ATT_REPEAT)
    new_seq_list[attr_pos - 1 : attr_pos - 1 + len(ATT_REPEAT)] = list(ATT_REPEAT)
    new_seq = "".join(new_seq_list)

    new_features: list[GBFeature] = []
    # Add explicit tRNA feature for trnS-Lys so att hotspot detection finds it cleanly
    new_features.append(
        GBFeature(
            feature_type="tRNA",
            location="87..162",
            start=87,
            end=162,
            strand="+",
            parts=[(87, 162)],
            qualifiers={"gene": ["trnS-Lys"], "product": ["tRNA-Lys"]},
        )
    )

    for feat in rec.features:
        if feat.feature_type == "source":
            continue
        s, e = feat.start, feat.end
        if s >= insert_pos:
            s += delta
            e += delta
        elif e >= insert_pos:
            e += delta
        new_features.append(
            GBFeature(
                feature_type=feat.feature_type,
                location=f"{s}..{e}",
                start=s,
                end=e,
                strand=feat.strand,
                parts=[(s, e)],
                qualifiers={k: list(v) for k, v in feat.qualifiers.items()},
            )
        )

    # Add inserted cargo CDS features inside ICEBs1
    cur = insert_pos + len(spacer)
    cargo_specs = [
        (
            "CARGO_TETM",
            "tet(M)",
            "tetracycline resistance ribosomal protection protein Tet(M)",
            TETM_PROT,
            tetm_dna,
        ),
        ("CARGO_CADA", "cadA", "cadmium-translocating P-type ATPase CadA", CADA_PROT, cada_dna),
        (
            "CARGO_RM_M",
            "hsdM",
            "Type II restriction-modification methyltransferase",
            RM_MT_PROT,
            rm_m_dna,
        ),
        ("CARGO_RM_R", "hsdR", "Type II restriction endonuclease", RM_RE_PROT, rm_r_dna),
    ]
    for lt, gene, product, prot, dna in cargo_specs:
        s = cur
        e = cur + len(dna) - 1
        new_features.append(
            GBFeature(
                feature_type="CDS",
                location=f"{s}..{e}",
                start=s,
                end=e,
                strand="+",
                parts=[(s, e)],
                qualifiers={
                    "locus_tag": [lt],
                    "gene": [gene],
                    "product": [product],
                    "protein_id": [f"WP_{lt}.1"],
                    "codon_start": ["1"],
                    "transl_table": ["11"],
                    "translation": [prot],
                },
            )
        )
        cur = e + 1 + len(spacer)

    new_features.sort(key=lambda f: (f.start, f.end))
    out_path = fixtures_dir / "mini_genome_cargo.gbff"
    gbk_text = format_genbank_record(
        locus_name="NC_000964_CARGO",
        definition=(
            "Bacillus subtilis ICEBs1 synthetic cargo test genome with tet(M), cadA, RM, attL/attR"
        ),
        organism="Bacillus subtilis subsp. subtilis str. 168",
        sequence=new_seq,
        features=new_features,
        topology="linear",
    )
    out_path.write_text(gbk_text, encoding="utf-8")
    return out_path


def build_mini_amrfinder_db(dbs_dir: Path, fixtures_dir: Path) -> Path:
    """Create minimal AMRFinderPlus v4 fixture DB covering tet(M), erm(B), cadA, and qacA."""
    db_dir = dbs_dir / "mini_amrfinder_db"
    db_dir.mkdir(parents=True, exist_ok=True)
    (db_dir / "version.txt").write_text("2026-01-01.1\n", encoding="utf-8")
    (db_dir / "database_format_version.txt").write_text("4.0\n", encoding="utf-8")

    entries = [
        (
            "WP_000691736.1|tet(M)|tet(M)|AMR|AMR|TETRACYCLINE|TETRACYCLINE|"
            "tetracycline resistance ribosomal protection protein Tet(M)",
            TETM_PROT,
        ),
        (
            "WP_001038795.1|erm(B)|erm(B)|AMR|AMR|MACROLIDE|MACROLIDE|"
            "23S rRNA (adenine(2058)-N(6))-methyltransferase Erm(B)",
            ERMB_PROT,
        ),
        (
            "WP_002354891.1|cadA|cadA|STRESS|METAL|CADMIUM|CADMIUM|"
            "cadmium-translocating P-type ATPase CadA",
            CADA_PROT,
        ),
        (
            "WP_000018312.1|qacA|qacA|STRESS|BIOCIDE|QUATERNARY_AMMONIUM|"
            "QUATERNARY_AMMONIUM|multidrug efflux MFS transporter QacA",
            QACA_PROT,
        ),
    ]
    faa_lines = [f">{hdr}\n{wrap_fasta(seq)}\n" for hdr, seq in entries]
    (db_dir / "AMRProt.fa").write_text("".join(faa_lines), encoding="utf-8")

    dna_entries = [
        ("WP_000691736.1|tet(M)|AMR|TETRACYCLINE", _back_translate(TETM_PROT)),
        ("WP_001038795.1|erm(B)|AMR|MACROLIDE", _back_translate(ERMB_PROT)),
        ("WP_002354891.1|cadA|STRESS|CADMIUM", _back_translate(CADA_PROT)),
        ("WP_000018312.1|qacA|STRESS|QUATERNARY_AMMONIUM", _back_translate(QACA_PROT)),
    ]
    dna_lines = [f">{hdr}\n{wrap_fasta(seq)}\n" for hdr, seq in dna_entries]
    (db_dir / "AMR_DNA-General.fa").write_text("".join(dna_lines), encoding="utf-8")

    template_hmm = fixtures_dir / "IntTyr_Phage_integrase.hmm"
    amr_lib_text = _build_valid_hmm(
        template_hmm,
        name="tet_M_fam",
        acc="NF033146.1",
        desc="tetracycline resistance ribosomal protection protein Tet(M)",
    ) + _build_valid_hmm(
        template_hmm,
        name="cadA_fam",
        acc="NF033201.1",
        desc="cadmium-translocating P-type ATPase CadA",
    )
    (db_dir / "AMR.LIB").write_text(amr_lib_text, encoding="utf-8")

    fam_rows = [
        "\t".join(
            [
                "fam_id",
                "parent_id",
                "gene_symbol",
                "hmm_id",
                "length",
                "complete_cutoff",
                "partial_cutoff",
                "scope",
                "element_type",
                "element_subtype",
                "class",
                "subclass",
                "name",
            ]
        ),
        "\t".join(["ALL", "", "-", "-", "0", "0", "0", "core", "AMR", "AMR", "-", "-", "root"]),
        "\t".join(
            [
                "AMR",
                "ALL",
                "-",
                "-",
                "0",
                "0",
                "0",
                "core",
                "AMR",
                "AMR",
                "-",
                "-",
                "Antimicrobial resistance",
            ]
        ),
        "\t".join(
            [
                "STRESS",
                "ALL",
                "-",
                "-",
                "0",
                "0",
                "0",
                "plus",
                "STRESS",
                "STRESS",
                "-",
                "-",
                "Stress response",
            ]
        ),
        "\t".join(
            [
                "TET",
                "AMR",
                "-",
                "-",
                "0",
                "0",
                "0",
                "core",
                "AMR",
                "AMR",
                "TETRACYCLINE",
                "TETRACYCLINE",
                "Tetracycline resistance",
            ]
        ),
        "\t".join(
            [
                "ERM",
                "AMR",
                "-",
                "-",
                "0",
                "0",
                "0",
                "core",
                "AMR",
                "AMR",
                "MACROLIDE",
                "MACROLIDE",
                "Macrolide resistance",
            ]
        ),
        "\t".join(
            [
                "METAL",
                "STRESS",
                "-",
                "-",
                "0",
                "0",
                "0",
                "plus",
                "STRESS",
                "METAL",
                "METAL",
                "METAL",
                "Heavy metal resistance",
            ]
        ),
        "\t".join(
            [
                "BIOCIDE",
                "STRESS",
                "-",
                "-",
                "0",
                "0",
                "0",
                "plus",
                "STRESS",
                "BIOCIDE",
                "BIOCIDE",
                "BIOCIDE",
                "Biocide resistance",
            ]
        ),
        "\t".join(
            [
                "tet(M)",
                "TET",
                "tet(M)",
                "NF033146.1",
                "240",
                "27.1",
                "27.1",
                "core",
                "AMR",
                "AMR",
                "TETRACYCLINE",
                "TETRACYCLINE",
                "tetracycline resistance ribosomal protection protein Tet(M)",
            ]
        ),
        "\t".join(
            [
                "erm(B)",
                "ERM",
                "erm(B)",
                "-",
                "160",
                "0",
                "0",
                "core",
                "AMR",
                "AMR",
                "MACROLIDE",
                "MACROLIDE",
                "23S rRNA methyltransferase Erm(B)",
            ]
        ),
        "\t".join(
            [
                "cadA",
                "METAL",
                "cadA",
                "NF033201.1",
                "230",
                "27.1",
                "27.1",
                "plus",
                "STRESS",
                "METAL",
                "CADMIUM",
                "CADMIUM",
                "cadmium-translocating P-type ATPase CadA",
            ]
        ),
        "\t".join(
            [
                "qacA",
                "BIOCIDE",
                "qacA",
                "-",
                "160",
                "0",
                "0",
                "plus",
                "STRESS",
                "BIOCIDE",
                "QUATERNARY_AMMONIUM",
                "QUATERNARY_AMMONIUM",
                "multidrug efflux MFS transporter QacA",
            ]
        ),
    ]
    (db_dir / "fam.tab").write_text("\n".join(fam_rows) + "\n", encoding="utf-8")

    hierarchy_rows = [
        "\t".join(
            [
                "node_id",
                "parent_node_id",
                "gene_symbol",
                "hmm_id",
                "hmm_tc1",
                "hmm_tc2",
                "blast_ident",
                "blast_cov",
                "scope",
                "type",
                "subtype",
                "class",
                "subclass",
                "name",
            ]
        ),
        "\t".join(
            [
                "ALL",
                "",
                "",
                "",
                "",
                "",
                "",
                "",
                "core",
                "AMR",
                "AMR",
                "",
                "",
                "Reference gene root",
            ]
        ),
        "\t".join(
            [
                "AMR",
                "ALL",
                "",
                "",
                "",
                "",
                "",
                "",
                "core",
                "AMR",
                "AMR",
                "",
                "",
                "Antimicrobial resistance",
            ]
        ),
        "\t".join(
            [
                "STRESS",
                "ALL",
                "",
                "",
                "",
                "",
                "",
                "",
                "plus",
                "STRESS",
                "STRESS",
                "",
                "",
                "Stress response",
            ]
        ),
        "\t".join(
            [
                "TET",
                "AMR",
                "",
                "",
                "",
                "",
                "",
                "",
                "core",
                "AMR",
                "AMR",
                "TETRACYCLINE",
                "TETRACYCLINE",
                "Tetracycline resistance",
            ]
        ),
        "\t".join(
            [
                "ERM",
                "AMR",
                "",
                "",
                "",
                "",
                "",
                "",
                "core",
                "AMR",
                "AMR",
                "MACROLIDE",
                "MACROLIDE",
                "Macrolide resistance",
            ]
        ),
        "\t".join(
            [
                "METAL",
                "STRESS",
                "",
                "",
                "",
                "",
                "",
                "",
                "plus",
                "STRESS",
                "METAL",
                "METAL",
                "METAL",
                "Heavy metal resistance",
            ]
        ),
        "\t".join(
            [
                "BIOCIDE",
                "STRESS",
                "",
                "",
                "",
                "",
                "",
                "",
                "plus",
                "STRESS",
                "BIOCIDE",
                "BIOCIDE",
                "BIOCIDE",
                "Biocide resistance",
            ]
        ),
        "\t".join(
            [
                "tet(M)",
                "TET",
                "tet(M)",
                "NF033146.1",
                "27.1",
                "27.1",
                "90",
                "90",
                "core",
                "AMR",
                "AMR",
                "TETRACYCLINE",
                "TETRACYCLINE",
                "tetracycline resistance ribosomal protection protein Tet(M)",
            ]
        ),
        "\t".join(
            [
                "erm(B)",
                "ERM",
                "erm(B)",
                "",
                "",
                "",
                "90",
                "90",
                "core",
                "AMR",
                "AMR",
                "MACROLIDE",
                "MACROLIDE",
                "23S rRNA (adenine(2058)-N(6))-methyltransferase Erm(B)",
            ]
        ),
        "\t".join(
            [
                "cadA",
                "METAL",
                "cadA",
                "NF033201.1",
                "27.1",
                "27.1",
                "90",
                "90",
                "plus",
                "STRESS",
                "METAL",
                "CADMIUM",
                "CADMIUM",
                "cadmium-translocating P-type ATPase CadA",
            ]
        ),
        "\t".join(
            [
                "qacA",
                "BIOCIDE",
                "qacA",
                "",
                "",
                "",
                "90",
                "90",
                "plus",
                "STRESS",
                "BIOCIDE",
                "QUATERNARY_AMMONIUM",
                "QUATERNARY_AMMONIUM",
                "multidrug efflux MFS transporter QacA",
            ]
        ),
    ]
    (db_dir / "ReferenceGeneHierarchy.txt").write_text(
        "\n".join(hierarchy_rows) + "\n", encoding="utf-8"
    )

    taxgroup_rows = [
        "#taxgroup\tgpipe_taxgroup\tnumber_of_nucl_ref_genes",
        "Bacillus_subtilis\tBacillus subtilis\t0",
        "Enterococcus_faecalis\tEnterococcus faecalis\t0",
        "Staphylococcus_aureus\tStaphylococcus aureus\t0",
    ]
    (db_dir / "taxgroup.tab").write_text("\n".join(taxgroup_rows) + "\n", encoding="utf-8")
    return db_dir


def build_mini_card_db(dbs_dir: Path) -> Path:
    """Create CARD v3.x fixture (ARO:3000186 tet(M), ARO:3000375 ErmB, ARO:3002639 vanA)."""
    db_dir = dbs_dir / "mini_card_db"
    db_dir.mkdir(parents=True, exist_ok=True)
    card_data = {
        "_version": "3.3.0-mini",
        "_timestamp": "2026-01-01T00:00:00Z",
        "_comment": "Miniature CARD v3.x fixture for offline advena testing",
        "1": {
            "model_id": "1",
            "model_name": "tet(M)",
            "model_type": "protein homolog model",
            "model_type_id": "40292",
            "ARO_accession": "3000186",
            "ARO_id": "36236",
            "ARO_name": "tet(M)",
            "ARO_description": (
                "Tet(M) is a ribosomal protection protein that confers tetracycline resistance."
            ),
            "cvterm_id": "36236",
            "ARO_category": {
                "3000050": {
                    "category_aro_accession": "3000050",
                    "category_aro_cvterm_id": "36100",
                    "category_aro_name": "tetracycline antibiotic",
                    "category_aro_class_name": "Drug Class",
                    "category_aro_description": "Tetracycline class of antibiotics.",
                },
                "3000172": {
                    "category_aro_accession": "3000172",
                    "category_aro_cvterm_id": "36222",
                    "category_aro_name": "antibiotic target protection",
                    "category_aro_class_name": "Resistance Mechanism",
                    "category_aro_description": "Protection of the antibiotic target site.",
                },
            },
            "model_sequences": {
                "sequence": {
                    "101": {
                        "protein_sequence": {"sequence": TETM_PROT},
                        "dna_sequence": {"sequence": _back_translate(TETM_PROT)},
                        "NCBI_taxonomy": {
                            "NCBI_taxonomy_id": "1423",
                            "NCBI_taxonomy_name": "Bacillus subtilis",
                        },
                    }
                }
            },
            "model_param": {
                "blastp_bit_score": {
                    "param_type": "BLASTP bit_score",
                    "param_type_id": "40310",
                    "param_value": "200",
                }
            },
        },
        "2": {
            "model_id": "2",
            "model_name": "ErmB",
            "model_type": "protein homolog model",
            "model_type_id": "40292",
            "ARO_accession": "3000375",
            "ARO_id": "36425",
            "ARO_name": "ErmB",
            "ARO_description": (
                "ErmB is a 23S rRNA methyltransferase conferring macrolide resistance."
            ),
            "cvterm_id": "36425",
            "ARO_category": {
                "3000002": {
                    "category_aro_accession": "3000002",
                    "category_aro_cvterm_id": "36052",
                    "category_aro_name": "macrolide antibiotic",
                    "category_aro_class_name": "Drug Class",
                    "category_aro_description": "Macrolide class of antibiotics.",
                },
                "3000175": {
                    "category_aro_accession": "3000175",
                    "category_aro_cvterm_id": "36225",
                    "category_aro_name": "antibiotic target alteration",
                    "category_aro_class_name": "Resistance Mechanism",
                    "category_aro_description": "Modification of the antibiotic target site.",
                },
            },
            "model_sequences": {
                "sequence": {
                    "102": {
                        "protein_sequence": {"sequence": ERMB_PROT},
                        "dna_sequence": {"sequence": _back_translate(ERMB_PROT)},
                        "NCBI_taxonomy": {
                            "NCBI_taxonomy_id": "1351",
                            "NCBI_taxonomy_name": "Enterococcus faecalis",
                        },
                    }
                }
            },
            "model_param": {
                "blastp_bit_score": {
                    "param_type": "BLASTP bit_score",
                    "param_type_id": "40310",
                    "param_value": "200",
                }
            },
        },
        "3": {
            "model_id": "3",
            "model_name": "vanA",
            "model_type": "protein homolog model",
            "model_type_id": "40292",
            "ARO_accession": "3002639",
            "ARO_id": "38689",
            "ARO_name": "vanA",
            "ARO_description": ("VanA is a D-Ala-D-Lac ligase conferring glycopeptide resistance."),
            "cvterm_id": "38689",
            "ARO_category": {
                "3000009": {
                    "category_aro_accession": "3000009",
                    "category_aro_cvterm_id": "36059",
                    "category_aro_name": "glycopeptide antibiotic",
                    "category_aro_class_name": "Drug Class",
                    "category_aro_description": "Glycopeptide class of antibiotics.",
                },
                "3000175": {
                    "category_aro_accession": "3000175",
                    "category_aro_cvterm_id": "36225",
                    "category_aro_name": "antibiotic target alteration",
                    "category_aro_class_name": "Resistance Mechanism",
                    "category_aro_description": "Modification of the antibiotic target site.",
                },
            },
            "model_sequences": {
                "sequence": {
                    "103": {
                        "protein_sequence": {"sequence": VANA_PROT},
                        "dna_sequence": {"sequence": _back_translate(VANA_PROT)},
                        "NCBI_taxonomy": {
                            "NCBI_taxonomy_id": "1352",
                            "NCBI_taxonomy_name": "Enterococcus faecium",
                        },
                    }
                }
            },
            "model_param": {
                "blastp_bit_score": {
                    "param_type": "BLASTP bit_score",
                    "param_type_id": "40310",
                    "param_value": "200",
                }
            },
        },
    }
    (db_dir / "card.json").write_text(json.dumps(card_data, indent=2) + "\n", encoding="utf-8")
    return db_dir


def build_mini_eggnog_db(dbs_dir: Path) -> Path:
    """Create eggNOG-mapper v2.1.12 compatible SQLite DBs (eggnog.db, eggnog.taxa.db) and FAA."""
    db_dir = dbs_dir / "mini_eggnog_db"
    db_dir.mkdir(parents=True, exist_ok=True)

    eggnog_db = db_dir / "eggnog.db"
    if eggnog_db.exists():
        eggnog_db.unlink()
    conn = sqlite3.connect(eggnog_db)
    conn.execute(
        """
        CREATE TABLE version (
            version VARCHAR(16) PRIMARY KEY
        ) WITHOUT ROWID
        """
    )
    conn.execute("INSERT INTO version VALUES (?)", ("5.0.2",))
    conn.execute(
        """
        CREATE TABLE og (
            og VARCHAR(16),
            level VARCHAR(16),
            nm INTEGER,
            description TEXT,
            COG_categories VARCHAR(8),
            PRIMARY KEY (og, level)
        ) WITHOUT ROWID
        """
    )
    conn.executemany(
        "INSERT INTO og VALUES (?, ?, ?, ?, ?)",
        [
            ("COG0582", "1", 500, "Site-specific recombinase / Tyrosine integrase", "L"),
            ("COG0582", "2", 450, "Site-specific recombinase / Tyrosine integrase", "L"),
            ("2IGK4", "1239", 120, "ICEBs1/Tn916 family tyrosine integrase", "L"),
            (
                "COG0480",
                "1",
                800,
                "Translation elongation factor / Tet(M) ribosomal protection",
                "J",
            ),
            (
                "COG0480",
                "2",
                700,
                "Translation elongation factor / Tet(M) ribosomal protection",
                "J",
            ),
            ("2IQ7R", "1239", 140, "Tetracycline resistance protein Tet(M)", "J"),
            ("COG2217", "1", 600, "Cation-transporting P-type ATPase CadA", "P"),
            ("COG2217", "2", 520, "Cation-transporting P-type ATPase CadA", "P"),
            ("2IT4N", "1239", 110, "Heavy metal translocating P-type ATPase", "P"),
            ("COG0863", "1", 400, "DNA modification methylase", "L"),
            ("COG0863", "2", 350, "DNA modification methylase", "L"),
            ("2IR9K", "1239", 95, "Type II restriction-modification system", "V"),
        ],
    )
    conn.execute(
        """
        CREATE TABLE event (
            i INTEGER PRIMARY KEY,
            level VARCHAR(16),
            og VARCHAR(16),
            side1 TEXT,
            side2 TEXT
        )
        """
    )
    conn.executemany(
        "INSERT INTO event VALUES (?, ?, ?, ?, ?)",
        [
            (1, "1239", "2IQ7R", "1423.CARGO_TETM", ""),
            (2, "1239", "2IT4N", "1423.CARGO_CADA", ""),
            (3, "1239", "2IR9K", "1423.CARGO_RM_M", ""),
        ],
    )
    conn.execute(
        """
        CREATE TABLE prots (
            name VARCHAR(32) PRIMARY KEY,
            bigg_reaction VARCHAR(32),
            gos TEXT,
            pfam TEXT,
            pname VARCHAR(32),
            ogs VARCHAR(256),
            orthoindex VARCHAR(256),
            kegg_ko VARCHAR(256),
            kegg_cog VARCHAR(256),
            kegg_disease VARCHAR(256),
            kegg_ec VARCHAR(256),
            kegg_brite VARCHAR(256),
            kegg_rclass VARCHAR(256),
            kegg_tc VARCHAR(256),
            kegg_cazy VARCHAR(256),
            kegg_pathway VARCHAR(256)
        ) WITHOUT ROWID
        """
    )
    conn.executemany(
        "INSERT INTO prots VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
        [
            (
                "1423.CARGO_TETM",
                "",
                "GO:0046677",
                "GTP_EFTU,Tet_M",
                "tetM",
                "COG0480@1|root,COG0480@2|Bacteria,2IQ7R@1239|Firmicutes",
                "1",
                "ko:K18220",
                "J",
                "",
                "",
                "ko00000",
                "",
                "",
                "",
                "ko01503",
            ),
            (
                "1423.CARGO_CADA",
                "",
                "GO:0046686",
                "E1-E2_ATPase,HMA",
                "cadA",
                "COG2217@1|root,COG2217@2|Bacteria,2IT4N@1239|Firmicutes",
                "2",
                "ko:K01534",
                "P",
                "",
                "7.2.2.21",
                "ko00000",
                "",
                "3.A.3.2",
                "",
                "ko02000",
            ),
            (
                "1423.CARGO_RM_M",
                "",
                "GO:0009307",
                "N6_Mtase",
                "hsdM",
                "COG0863@1|root,COG0863@2|Bacteria,2IR9K@1239|Firmicutes",
                "3",
                "ko:K03427",
                "V",
                "",
                "2.1.1.72",
                "ko00000",
                "",
                "",
                "",
                "",
            ),
        ],
    )
    conn.commit()
    conn.close()

    taxa_db = db_dir / "eggnog.taxa.db"
    if taxa_db.exists():
        taxa_db.unlink()
    tconn = sqlite3.connect(taxa_db)
    tconn.execute("CREATE TABLE stats (version INT PRIMARY KEY)")
    tconn.execute("INSERT INTO stats VALUES (2)")
    tconn.execute(
        """
        CREATE TABLE species (
            taxid INT PRIMARY KEY,
            parent INT,
            spname VARCHAR,
            common VARCHAR,
            rank VARCHAR,
            track TEXT
        )
        """
    )
    tconn.executemany(
        "INSERT INTO species VALUES (?, ?, ?, ?, ?, ?)",
        [
            (1, 1, "root", "", "no rank", "1"),
            (2, 1, "Bacteria", "", "superkingdom", "2,1"),
            (1239, 2, "Firmicutes", "", "phylum", "1239,2,1"),
            (91061, 1239, "Bacilli", "", "class", "91061,1239,2,1"),
            (1385, 91061, "Bacillales", "", "order", "1385,91061,1239,2,1"),
            (186817, 1385, "Bacillaceae", "", "family", "186817,1385,91061,1239,2,1"),
            (1386, 186817, "Bacillus", "", "genus", "1386,186817,1385,91061,1239,2,1"),
            (
                1423,
                1386,
                "Bacillus subtilis",
                "",
                "species",
                "1423,1386,186817,1385,91061,1239,2,1",
            ),
        ],
    )
    tconn.execute(
        """
        CREATE TABLE synonym (
            taxid INT,
            spname VARCHAR,
            PRIMARY KEY (spname, taxid)
        )
        """
    )
    tconn.executemany(
        "INSERT INTO synonym VALUES (?, ?)",
        [
            (1239, "Bacillota"),
            (1423, "Bacillus subtilis subsp. subtilis"),
        ],
    )
    tconn.execute(
        """
        CREATE TABLE merged (
            taxid_old INT PRIMARY KEY,
            taxid_new INT
        )
        """
    )
    tconn.execute("INSERT INTO merged VALUES (?, ?)", (224308, 1423))
    tconn.commit()
    tconn.close()

    (db_dir / "eggnog_proteins.faa").write_text(
        f">1423.CARGO_TETM\n{wrap_fasta(TETM_PROT)}\n"
        f">1423.CARGO_CADA\n{wrap_fasta(CADA_PROT)}\n"
        f">1423.CARGO_RM_M\n{wrap_fasta(RM_MT_PROT)}\n",
        encoding="utf-8",
    )
    return db_dir


def build_mini_defensefinder_db(dbs_dir: Path, fixtures_dir: Path) -> Path:
    """Create minimal DefenseFinder/MacSyFinder model tree with valid HMMER3/f profiles."""
    db_dir = dbs_dir / "mini_defensefinder_db"
    defs_dir = db_dir / "defense-finder-models" / "definitions" / "RM"
    prof_dir = db_dir / "defense-finder-models" / "profiles"
    defs_dir.mkdir(parents=True, exist_ok=True)
    prof_dir.mkdir(parents=True, exist_ok=True)

    (db_dir / "defense-finder-models" / "metadata.yml").write_text(
        "vers: '2.0.0-mini'\n"
        "format_version: 2\n"
        "short_desc: 'Miniature offline DefenseFinder models for advena testing'\n"
        "maintainer:\n"
        "  name: 'advena-test'\n"
        "  email: 'advena-test@example.org'\n"
        "cite:\n"
        "  - 'Tesson et al., Nature Communications 2022'\n"
        "doc: 'https://github.com/mdmparis/defense-finder-models'\n"
        "license: 'GPL-3.0'\n"
        "copyright: 'advena-test'\n",
        encoding="utf-8",
    )
    (defs_dir / "RM_Type_II.xml").write_text(
        '<model inter_gene_max_space="5" min_mandatory_genes_required="1" '
        'min_genes_required="2" vers="2.0">\n'
        '  <gene name="RM_Type_II__Type_II_MTases" presence="mandatory"/>\n'
        '  <gene name="RM_Type_II__Type_II_REases" presence="accessory"/>\n'
        "</model>\n",
        encoding="utf-8",
    )
    template_hmm = fixtures_dir / "IntTyr_Phage_integrase.hmm"
    (prof_dir / "RM_Type_II__Type_II_MTases.hmm").write_text(
        _build_valid_hmm(
            template_hmm,
            name="RM_Type_II__Type_II_MTases",
            desc="Type II restriction-modification methyltransferase (mini fixture)",
        ),
        encoding="utf-8",
    )
    (prof_dir / "RM_Type_II__Type_II_REases.hmm").write_text(
        _build_valid_hmm(
            template_hmm,
            name="RM_Type_II__Type_II_REases",
            desc="Type II restriction-modification endonuclease (mini fixture)",
        ),
        encoding="utf-8",
    )
    return db_dir


def build_all(fixtures_dir: Path) -> None:
    dbs_dir = fixtures_dir / "dbs"
    dbs_dir.mkdir(parents=True, exist_ok=True)
    build_mini_genome_cargo(fixtures_dir)
    build_mini_amrfinder_db(dbs_dir, fixtures_dir)
    build_mini_card_db(dbs_dir)
    build_mini_eggnog_db(dbs_dir)
    build_mini_defensefinder_db(dbs_dir, fixtures_dir)

    samplesheet = fixtures_dir / "samplesheet_full.csv"
    samplesheet.write_text(
        "sample,genbank\nmini_cargo,mini_genome_cargo.gbff\nmini_ref,mini_genome.gbff\n",
        encoding="utf-8",
    )


def main() -> int:
    fixtures_dir = REPO_ROOT / "tests" / "fixtures"
    build_all(fixtures_dir)
    print(f"Built miniature fixture databases and cargo genome under {fixtures_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
