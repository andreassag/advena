#!/usr/bin/env python3
"""Validate advena JSON schemas against Draft 2020-12 and check nextflow.config parity."""

from __future__ import annotations

import json
import re
import sys
from pathlib import Path

import jsonschema

REPO_ROOT = Path(__file__).resolve().parents[2]


def extract_config_params(config_path: Path) -> set[str]:
    """Extract top-level parameter names from the params { ... } block in nextflow.config."""
    text = config_path.read_text(encoding="utf-8")
    start_match = re.search(r"^\s*params\s*\{", text, re.MULTILINE)
    if not start_match:
        raise ValueError(f"Could not locate params block in {config_path}")

    depth = 0
    block_lines: list[str] = []
    started = False
    for line in text[start_match.start() :].splitlines():
        depth += line.count("{") - line.count("}")
        if started and depth >= 1:
            block_lines.append(line)
        if not started and "{" in line:
            started = True
        if started and depth == 0:
            break

    params: set[str] = set()
    for line in block_lines:
        stripped = re.sub(r"//.*$", "", line).strip()
        if not stripped:
            continue
        match = re.match(r"^([a-zA-Z_][a-zA-Z0-9_]*)\s*=", stripped)
        if match:
            params.add(match.group(1))
    return params


def extract_schema_params(schema: dict) -> set[str]:
    """Collect all parameter names from top-level properties and $defs."""
    params: set[str] = set(schema.get("properties", {}).keys())
    for definition in schema.get("$defs", {}).values():
        if isinstance(definition, dict):
            params.update(definition.get("properties", {}).keys())
    return params


def validate_all(repo_root: Path = REPO_ROOT) -> list[str]:
    """Validate schemas and return a list of error messages (empty on success)."""
    errors: list[str] = []

    schema_paths = [
        repo_root / "nextflow_schema.json",
        repo_root / "assets" / "schema_input.json",
    ]
    loaded: dict[str, dict] = {}

    for path in schema_paths:
        if not path.exists():
            errors.append(f"Missing schema file: {path}")
            continue
        try:
            schema = json.loads(path.read_text(encoding="utf-8"))
        except json.JSONDecodeError as exc:
            errors.append(f"Invalid JSON in {path}: {exc}")
            continue
        try:
            jsonschema.Draft202012Validator.check_schema(schema)
        except jsonschema.SchemaError as exc:
            errors.append(f"Draft 2020-12 schema error in {path}: {exc.message}")
            continue
        loaded[path.name] = schema

    nf_schema = loaded.get("nextflow_schema.json")
    if nf_schema is not None:
        required = {"$schema", "$id", "title", "description", "type", "allOf"}
        missing = required - nf_schema.keys()
        if missing:
            errors.append(
                f"nextflow_schema.json missing required top-level keys: {sorted(missing)}"
            )
        if "properties" not in nf_schema and "$defs" not in nf_schema:
            errors.append("nextflow_schema.json must define 'properties' or '$defs'")

        defs = nf_schema.get("$defs", {})
        for entry in nf_schema.get("allOf", []):
            ref = entry.get("$ref", "")
            if ref.startswith("#/$defs/"):
                key = ref.split("#/$defs/", 1)[1]
                if key not in defs:
                    errors.append(f"Unresolved $ref in nextflow_schema.json: {ref}")

        config_path = repo_root / "nextflow.config"
        if config_path.exists():
            cfg_params = extract_config_params(config_path)
            sch_params = extract_schema_params(nf_schema)
            missing_in_schema = sorted(cfg_params - sch_params)
            missing_in_config = sorted(sch_params - cfg_params)
            if missing_in_schema:
                errors.append(
                    f"Parameters in nextflow.config missing from nextflow_schema.json: "
                    f"{missing_in_schema}"
                )
            if missing_in_config:
                errors.append(
                    f"Parameters in nextflow_schema.json missing from nextflow.config: "
                    f"{missing_in_config}"
                )

    input_schema = loaded.get("schema_input.json")
    if input_schema is not None:
        validator = jsonschema.Draft202012Validator(input_schema)
        valid_cases = [
            [{"sample": "s1", "genbank": "a.gbff"}],
            [{"sample": "s1", "genbank": "a.gbk.gz"}],
            [{"sample": "s2", "fasta": "b.fna"}],
            [{"sample": "s2", "fasta": "b.fasta.gz"}],
            [{"sample": "s3", "genbank": "a.gb", "fasta": "b.fa"}],
        ]
        for case in valid_cases:
            errs = list(validator.iter_errors(case))
            if errs:
                errors.append(f"schema_input.json unexpectedly rejected {case}: {errs[0].message}")

        invalid_cases = [
            [{"sample": "s1"}],
            [{"sample": "bad sample", "genbank": "a.gbff"}],
            [{"sample": "s1", "genbank": "a.txt"}],
            [{"sample": "s1", "fasta": "a.txt"}],
            [{"sample": "s1", "genbank": "a.gbff", "extra": "nope"}],
        ]
        for case in invalid_cases:
            if validator.is_valid(case):
                errors.append(f"schema_input.json unexpectedly accepted invalid case: {case}")

    return errors


def main() -> int:
    errors = validate_all()
    if errors:
        for err in errors:
            print(f"ERROR: {err}", file=sys.stderr)
        return 1
    print("All schemas are valid and synchronized with nextflow.config.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
