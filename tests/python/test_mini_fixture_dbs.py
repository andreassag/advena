"""Unit tests validating miniature offline fixture databases and static pipeline invariants."""

from __future__ import annotations

import json
import re
import sqlite3
import subprocess
import sys
import xml.etree.ElementTree as ET
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
DBS_DIR = REPO_ROOT / "tests" / "fixtures" / "dbs"


def _parse_hmmer3_profiles(hmm_path: Path) -> list[dict[str, object]]:
    """Parse and structurally validate one or more HMMER3/f profiles in a file."""
    text = hmm_path.read_text(encoding="utf-8")
    raw_blocks = [b.strip() for b in text.split("\n//") if b.strip()]
    assert raw_blocks, f"No HMMER3/f profile blocks found in {hmm_path}"
    profiles: list[dict[str, object]] = []

    for block in raw_blocks:
        lines = block.splitlines()
        assert lines[0].startswith("HMMER3/f"), f"Invalid HMMER3 header in {hmm_path}: {lines[0]}"
        meta: dict[str, str] = {}
        hmm_idx = -1
        for i, line in enumerate(lines):
            if line.startswith("HMM "):
                hmm_idx = i
                break
            parts = line.split(None, 1)
            if len(parts) == 2:
                meta[parts[0]] = parts[1].strip()

        assert "NAME" in meta, f"Missing NAME in {hmm_path}"
        assert "LENG" in meta and meta["LENG"].isdigit(), f"Missing/invalid LENG in {hmm_path}"
        assert meta.get("ALPH") == "amino", f"Expected ALPH amino in {hmm_path}"
        assert hmm_idx != -1, f"Missing HMM emission matrix header in {hmm_path}"

        leng = int(meta["LENG"])
        # After HMM header line: transition label line, COMPO + 2 lines, then 3 lines per match node
        body_lines = lines[hmm_idx + 1 :]
        assert len(body_lines) >= 4 + 3 * leng, (
            f"Truncated HMM matrix in {hmm_path}: expected >= {4 + 3 * leng} lines, "
            f"found {len(body_lines)}"
        )
        assert body_lines[1].strip().startswith("COMPO"), f"Missing COMPO line in {hmm_path}"
        # Verify match node 1 and final match node LENG exist
        first_node = body_lines[4].strip().split()[0]
        last_node = body_lines[4 + 3 * (leng - 1)].strip().split()[0]
        assert first_node == "1", f"Expected first match node 1, got {first_node}"
        assert last_node == str(leng), f"Expected last match node {leng}, got {last_node}"
        profiles.append({"name": meta["NAME"], "acc": meta.get("ACC", ""), "leng": leng})

    return profiles


def test_mini_defensefinder_db_structure() -> None:
    root = DBS_DIR / "mini_defensefinder_db" / "defense-finder-models"
    meta_path = root / "metadata.yml"
    assert meta_path.is_file()
    meta_text = meta_path.read_text(encoding="utf-8")
    for key in (
        "vers:",
        "format_version:",
        "short_desc:",
        "maintainer:",
        "  name:",
        "  email:",
        "cite:",
        "doc:",
        "license:",
    ):
        assert key in meta_text, f"Missing DefenseFinder metadata.yml key: {key}"

    xml_path = root / "definitions" / "RM" / "RM_Type_II.xml"
    assert xml_path.is_file()

    tree = ET.parse(xml_path)
    genes = [elem.attrib["name"] for elem in tree.getroot().findall("gene")]
    assert set(genes) == {"RM_Type_II__Type_II_MTases", "RM_Type_II__Type_II_REases"}

    for gene in genes:
        hmm_file = root / "profiles" / f"{gene}.hmm"
        assert hmm_file.is_file(), f"Missing DefenseFinder profile: {hmm_file}"
        profiles = _parse_hmmer3_profiles(hmm_file)
        assert len(profiles) == 1
        assert profiles[0]["name"] == gene


def test_mini_amrfinder_db_structure() -> None:
    db = DBS_DIR / "mini_amrfinder_db"
    for fname in (
        "version.txt",
        "database_format_version.txt",
        "AMRProt.fa",
        "AMR_DNA-General.fa",
        "AMR.LIB",
        "fam.tab",
        "ReferenceGeneHierarchy.txt",
        "taxgroup.tab",
    ):
        assert (db / fname).is_file(), f"Missing AMRFinderPlus DB file: {fname}"
        assert (db / fname).stat().st_size > 0, f"Empty AMRFinderPlus DB file: {fname}"

    assert (db / "database_format_version.txt").read_text(encoding="utf-8").strip() == "4.0"

    profiles = _parse_hmmer3_profiles(db / "AMR.LIB")
    assert len(profiles) == 2
    accs = {p["acc"] for p in profiles}
    assert accs == {"NF033146.1", "NF033201.1"}

    prot_text = (db / "AMRProt.fa").read_text(encoding="utf-8")
    assert ">WP_000691736.1|tet(M)" in prot_text
    assert ">WP_002354891.1|cadA" in prot_text

    dna_text = (db / "AMR_DNA-General.fa").read_text(encoding="utf-8")
    assert ">WP_000691736.1|tet(M)" in dna_text

    hier_lines = (db / "ReferenceGeneHierarchy.txt").read_text(encoding="utf-8").splitlines()
    assert hier_lines[0].startswith("node_id\tparent_node_id\tgene_symbol\thmm_id")
    assert any("\ttet(M)\tNF033146.1\t" in line for line in hier_lines[1:])

    tax_lines = (db / "taxgroup.tab").read_text(encoding="utf-8").splitlines()
    assert tax_lines[0].startswith("#taxgroup\tgpipe_taxgroup")
    assert any(line.startswith("Bacillus_subtilis\t") for line in tax_lines[1:])


def test_mini_card_db_v3_schema() -> None:
    card_path = DBS_DIR / "mini_card_db" / "card.json"
    assert card_path.is_file()
    data = json.loads(card_path.read_text(encoding="utf-8"))

    for key in ("_version", "_timestamp", "_comment"):
        assert key in data and isinstance(data[key], str) and data[key]

    models = {k: v for k, v in data.items() if not k.startswith("_")}
    assert set(models.keys()) == {"1", "2", "3"}

    required_model_keys = {
        "model_id",
        "model_name",
        "model_type",
        "model_type_id",
        "ARO_accession",
        "ARO_name",
        "ARO_description",
        "cvterm_id",
        "ARO_category",
        "model_sequences",
        "model_param",
    }
    for model_id, entry in models.items():
        missing = required_model_keys - set(entry.keys())
        assert not missing, f"Model {model_id} missing CARD keys: {missing}"
        assert entry["model_type_id"] == "40292"
        for cat in entry["ARO_category"].values():
            assert "category_aro_cvterm_id" in cat
            assert "category_aro_accession" in cat
            assert "category_aro_name" in cat
            assert "category_aro_class_name" in cat
        bitscore = entry["model_param"]["blastp_bit_score"]
        assert bitscore["param_type"] == "BLASTP bit_score"
        assert int(bitscore["param_value"]) > 0
        seqs = entry["model_sequences"]["sequence"]
        assert len(seqs) >= 1


def test_mini_eggnog_db_v2_schema() -> None:
    db_dir = DBS_DIR / "mini_eggnog_db"
    eggnog_db = db_dir / "eggnog.db"
    taxa_db = db_dir / "eggnog.taxa.db"
    faa = db_dir / "eggnog_proteins.faa"
    assert eggnog_db.is_file() and taxa_db.is_file() and faa.is_file()

    with sqlite3.connect(eggnog_db) as conn:
        tables = {
            row[0]
            for row in conn.execute("SELECT name FROM sqlite_master WHERE type='table'").fetchall()
        }
        assert {"version", "og", "event", "prots"}.issubset(tables)

        ver = conn.execute("SELECT version FROM version").fetchone()
        assert ver is not None and ver[0] == "5.0.2"

        og_cols = {row[1] for row in conn.execute("PRAGMA table_info(og)").fetchall()}
        assert {"og", "level", "nm", "description", "COG_categories"} == og_cols

        prot_cols = {row[1] for row in conn.execute("PRAGMA table_info(prots)").fetchall()}
        assert {"name", "ogs", "gos", "pfam", "pname", "kegg_ko"}.issubset(prot_cols)
        prots = {row[0] for row in conn.execute("SELECT name FROM prots").fetchall()}
        assert {"1423.CARGO_TETM", "1423.CARGO_CADA", "1423.CARGO_RM_M"} == prots

        events = conn.execute("SELECT count(*) FROM event").fetchone()[0]
        assert events >= 3

    with sqlite3.connect(taxa_db) as tconn:
        ttables = {
            row[0]
            for row in tconn.execute("SELECT name FROM sqlite_master WHERE type='table'").fetchall()
        }
        assert {"stats", "species", "synonym", "merged"}.issubset(ttables)

        sp_cols = {row[1] for row in tconn.execute("PRAGMA table_info(species)").fetchall()}
        assert {"taxid", "parent", "spname", "common", "rank", "track"} == sp_cols

        bsu = tconn.execute("SELECT spname, track FROM species WHERE taxid = 1423").fetchone()
        assert bsu is not None and bsu[0] == "Bacillus subtilis"
        assert "1239" in bsu[1].split(",") and bsu[1].endswith(",1")

        syn = tconn.execute("SELECT taxid FROM synonym WHERE spname = ?", ("Bacillota",)).fetchone()
        assert syn is not None and syn[0] == 1239


def test_static_nextflow_guards_and_eval_sanitization() -> None:
    """Verify eval() version directives, AMRFinderPlus/geNomad guards, and MultiQC versions."""
    modules_dir = REPO_ROOT / "modules" / "local"
    for nf_file in sorted(modules_dir.rglob("main.nf")):
        text = nf_file.read_text(encoding="utf-8")
        for lineno, line in enumerate(text.splitlines(), start=1):
            if "eval(" in line:
                assert "| head -" not in line, f"SIGPIPE-prone head in {nf_file}:{lineno}"
                if "2>&1" in line:
                    assert "sed -n" in line, f"Unsanitized 2>&1 eval in {nf_file}:{lineno}"

    # Verify MULTIQC does not publish to topic: versions (prevents topic cycle/deadlock)
    mqc_text = (modules_dir / "multiqc" / "main.nf").read_text(encoding="utf-8")
    assert "topic: versions" not in mqc_text
    assert "emit: versions_multiqc" in mqc_text

    # Verify workflows/advena.nf subscribes to channel.topic("versions") once, forks via multiMap,
    # wires software_versions_mqc_versions.yml into MULTIQC, and emits versions_blocks
    advena_nf = (REPO_ROOT / "workflows" / "advena.nf").read_text(encoding="utf-8")
    assert advena_nf.count('channel.topic("versions")') == 1
    assert ".multiMap {" in advena_nf
    assert "software_versions_mqc_versions.yml" in advena_nf
    assert "sanitizeVersionValue" in advena_nf
    assert "sanitizeYamlKey" in advena_nf
    assert "versions_multiqc  = COHORT_COMPARATIVE.out.versions_multiqc" in advena_nf
    assert "versions_blocks   = ch_versions_blocks" in advena_nf

    # Verify main.nf writes pipeline_info/software_versions.yml from ADVENA.out.versions_blocks
    # without double-subscribing to channel.topic("versions") or cross-file def imports
    main_nf = (REPO_ROOT / "main.nf").read_text(encoding="utf-8")
    assert 'topic("versions")' not in main_nf
    assert "ADVENA.out.versions_blocks" in main_nf
    assert 'name: "software_versions.yml"' in main_nf

    # Verify AMRFINDERPLUS_RUN mutually exclusive input-mode branches + GFF feature check + hmmpress
    amr_nf = (modules_dir / "amrfinderplus" / "run" / "main.nf").read_text(encoding="utf-8")
    assert "if [ -s query.faa ]; then" in amr_nf
    assert "grep -q '^[^#]' query.gff3" in amr_nf
    assert (
        'if [ "\\${HAS_FAA}" -eq 1 ] && [ "\\${HAS_FNA}" -eq 1 ] && [ "\\${HAS_GFF}" -eq 1 ]; then'
        in amr_nf
    )
    assert "hmmpress" in amr_nf
    assert re.search(r"--protein query\.faa\s*\\", amr_nf) is None

    # Verify GENOMAD_ENDTOEND safely resolves nested genomad_db directory via -type f
    genomad_nf = (modules_dir / "genomad" / "endtoend" / "main.nf").read_text(encoding="utf-8")
    assert '-type f -name "genomad_db"' in genomad_nf


def test_amrfinderplus_input_mode_guard_shell(tmp_path: Path) -> None:
    """Execute AMRFINDERPLUS_RUN's exact input-selection bash logic under set -euo pipefail."""
    full_args = (
        "--protein query.faa --nucleotide query.fna --gff query.gff3 --annotation_format standard"
    )
    script_body = f"""
set -euo pipefail
if [ ! -s query.faa ] && [ ! -s query.fna ]; then
    echo "SKIP"
    exit 0
fi
HAS_FAA=0
HAS_FNA=0
HAS_GFF=0
if [ -s query.faa ]; then
    HAS_FAA=1
fi
if [ -s query.fna ]; then
    HAS_FNA=1
fi
if [ -s query.gff3 ] && grep -q '^[^#]' query.gff3; then
    HAS_GFF=1
fi
if [ "${{HAS_FAA}}" -eq 1 ] && [ "${{HAS_FNA}}" -eq 1 ] && [ "${{HAS_GFF}}" -eq 1 ]; then
    INPUT_ARGS="{full_args}"
elif [ "${{HAS_FAA}}" -eq 1 ]; then
    INPUT_ARGS="--protein query.faa"
else
    INPUT_ARGS="--nucleotide query.fna"
fi
echo "${{INPUT_ARGS}}"
"""
    cases = [
        ("", "", "", "SKIP"),
        ("", ">e1\nACGT\n", "##gff-version 3\n", "--nucleotide query.fna"),
        (">p1\nMKK\n", "", "", "--protein query.faa"),
        # Non-empty faa + fna with header-only gff3 MUST NOT pass -p + -n without -g
        (">p1\nMKK\n", ">e1\nACGT\n", "##gff-version 3\n", "--protein query.faa"),
        (
            ">p1\nMKK\n",
            ">e1\nACGT\n",
            "##gff-version 3\ne1\tadvena\tCDS\t1\t9\t.\t+\t0\tID=p1\n",
            full_args,
        ),
    ]
    for idx, (faa, fna, gff, expected) in enumerate(cases):
        case_dir = tmp_path / f"case_{idx}"
        case_dir.mkdir()
        (case_dir / "query.faa").write_text(faa, encoding="utf-8")
        (case_dir / "query.fna").write_text(fna, encoding="utf-8")
        (case_dir / "query.gff3").write_text(gff, encoding="utf-8")
        res = subprocess.run(
            ["/bin/bash", "-c", script_body],
            cwd=case_dir,
            capture_output=True,
            text=True,
            check=True,
        )
        assert res.stdout.strip() == expected


def test_empty_cohort_multiqc_yaml_valid(tmp_path: Path) -> None:
    """Verify build_cohort_summary emits valid YAML mapping (data: {}) even with 0 samples."""
    sys.path.insert(0, str(REPO_ROOT / "bin"))
    from build_cohort_summary import build_cohort_outputs

    n_me, n_og = build_cohort_outputs([], [], [], None, None, tmp_path)
    assert (n_me, n_og) == (0, 0)
    yaml_text = (tmp_path / "icescreen_summary_mqc.yaml").read_text(encoding="utf-8")
    assert "data: {}\n" in yaml_text
    assert "\ndata:\n" not in yaml_text


def _extract_block(text: str, header: str) -> str:
    """Extract an indented Nextflow section block (input:, output:, take:, emit:)."""
    lines = text.splitlines()
    collecting = False
    collected: list[str] = []
    for line in lines:
        stripped = line.strip()
        if not collecting:
            if stripped == header:
                collecting = True
            continue
        if (
            re.match(r"^(input|output|when|script|stub|take|main|emit):$", stripped)
            or stripped == "}"
        ):
            break
        collected.append(line)
    return "\n".join(collected)


def test_all_local_modules_and_subworkflows_have_valid_nf_tests() -> None:
    """Static audit: 100% nf-test coverage and input/output consistency across all modules."""
    # 1. Local modules
    module_nfs = sorted((REPO_ROOT / "modules" / "local").rglob("main.nf"))
    assert len(module_nfs) >= 45
    for nf_path in module_nfs:
        test_path = nf_path.parent / "tests" / "main.nf.test"
        assert test_path.is_file(), f"Missing nf-test for module {nf_path.parent}"
        nf_text = nf_path.read_text(encoding="utf-8")
        test_text = test_path.read_text(encoding="utf-8")

        m_proc = re.search(r"^\s*process\s+(\w+)\s*\{", nf_text, re.MULTILINE)
        assert m_proc, f"Could not parse process name in {nf_path}"
        proc_name = m_proc.group(1)

        assert f'process "{proc_name}"' in test_text, (
            f"{test_path} does not match process {proc_name}"
        )
        assert 'script "../main.nf"' in test_text
        assert 'tag "stub"' in test_text
        assert test_text.count('"""') % 2 == 0, f"Unbalanced triple quotes in {test_path}"

        # Count top-level input channels declared at 4-space indentation in input: block
        input_block = _extract_block(nf_text, "input:")
        declared_inputs = [
            ln
            for ln in input_block.splitlines()
            if re.match(r"^    (tuple|path|val|env|stdin)\b", ln)
        ]
        # Check referenced input indices in the first test block
        first_test = test_text.split("test(")[1]
        referenced_indices = {int(i) for i in re.findall(r"\binput\[(\d+)\]\s*=", first_test)}
        if declared_inputs:
            assert referenced_indices == set(range(len(declared_inputs))), (
                f"Input count mismatch in {test_path}: module declares {len(declared_inputs)} "
                f"inputs, test sets {sorted(referenced_indices)}"
            )
        else:
            assert not referenced_indices, (
                f"Module {nf_path} has 0 inputs but {test_path} sets {sorted(referenced_indices)}"
            )

        # Verify every referenced process.out.<emit> exists in output: block
        output_block = _extract_block(nf_text, "output:")
        emits = set(re.findall(r"\bemit:\s*(\w+)", output_block))
        for ref_emit in re.findall(r"\bprocess\.out\.(\w+)", test_text):
            assert ref_emit in emits, (
                f"{test_path} references process.out.{ref_emit} not in {sorted(emits)}"
            )

    # 2. Local subworkflows
    sw_nfs = sorted((REPO_ROOT / "subworkflows" / "local").rglob("main.nf"))
    assert len(sw_nfs) == 6
    for nf_path in sw_nfs:
        test_path = nf_path.parent / "tests" / "main.nf.test"
        assert test_path.is_file(), f"Missing nf-test for subworkflow {nf_path.parent}"
        nf_text = nf_path.read_text(encoding="utf-8")
        test_text = test_path.read_text(encoding="utf-8")

        m_wf = re.search(r"^\s*workflow\s+(\w+)\s*\{", nf_text, re.MULTILINE)
        assert m_wf, f"Could not parse workflow name in {nf_path}"
        wf_name = m_wf.group(1)

        assert f'workflow "{wf_name}"' in test_text, (
            f"{test_path} does not match workflow {wf_name}"
        )
        assert 'script "../main.nf"' in test_text
        assert 'tag "stub"' in test_text
        assert test_text.count('"""') % 2 == 0, f"Unbalanced triple quotes in {test_path}"

        take_block = _extract_block(nf_text, "take:")
        take_vars = [
            ln.strip().split()[0]
            for ln in take_block.splitlines()
            if ln.strip() and not ln.strip().startswith("//")
        ]
        first_test = test_text.split("test(")[1]
        referenced_indices = {int(i) for i in re.findall(r"\binput\[(\d+)\]\s*=", first_test)}
        assert referenced_indices == set(range(len(take_vars))), (
            f"Subworkflow input count mismatch in {test_path}: declares {len(take_vars)} "
            f"take vars, test sets {sorted(referenced_indices)}"
        )

        emit_block = _extract_block(nf_text, "emit:")
        emits = set(re.findall(r"^\s*(\w+)\s*=", emit_block, re.MULTILINE))
        for ref_emit in re.findall(r"\bworkflow\.out\.(\w+)", test_text):
            assert ref_emit in emits, (
                f"{test_path} references workflow.out.{ref_emit} not in {sorted(emits)}"
            )

        for cfg in re.findall(r'^\s*config\s+"([^"]+)"', test_text, re.MULTILINE):
            resolved = (test_path.parent / cfg).resolve()
            assert resolved.is_file(), f"Missing config path {cfg} referenced in {test_path}"

    # 3. Pipeline tests
    for ptest in sorted((REPO_ROOT / "tests").glob("*.nf.test")):
        ptext = ptest.read_text(encoding="utf-8")
        assert ptext.count('"""') % 2 == 0, f"Unbalanced triple quotes in {ptest}"
        for cfg in re.findall(r'^\s*config\s+"([^"]+)"', ptext, re.MULTILINE):
            resolved = (ptest.parent / cfg).resolve()
            assert resolved.is_file(), f"Missing config path {cfg} referenced in {ptest}"

    # 4. Guard against Groovy GString quote/newline unescaping bugs inside """..."""
    for tf in sorted(REPO_ROOT.rglob("*.nf.test")):
        if "modules/nf-core" in tf.as_posix():
            continue
        ttext = tf.read_text(encoding="utf-8")
        for m in re.finditer(r'"""(.*?)"""', ttext, re.DOTALL):
            body = m.group(1)
            for line in body.splitlines():
                assert not re.search(r'(?<!\\)\\"', line), (
                    f"Single-escaped quote inside Groovy GString in {tf}: {line.strip()}"
                )
                assert not re.search(r"(?<!\\)\\[ntr]", line), (
                    f"Single-escaped newline/tab inside Groovy GString in {tf}: {line.strip()}"
                )

    run_all_text = (REPO_ROOT / "tests" / "run_all.nf.test").read_text(encoding="utf-8")
    assert 'cargo_scope = "genome"' in run_all_text
    assert "strict_qc = true" in run_all_text
    assert "filter_non_bacillota = true" in run_all_text
    assert 'config "../conf/test_full.config"' in run_all_text
    assert 'profile "test"' not in run_all_text
