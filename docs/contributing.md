# Contributing

Contributions are welcome. This document describes how to set up a development environment and run the 4-tier local
verification harness before pushing changes.

---

## Development setup

```bash
git clone https://github.com/andreassag/advena.git
cd advena

# Install Node.js dependencies (Prettier)
npm install

# Install Python dev dependencies (pytest, ruff, jsonschema, pre-commit)
pip install -e ".[dev]"

# Install gh-act extension for local GitHub Actions verification
gh extension install https://github.com/nektos/gh-act

# Install pre-commit and git pre-push hooks
make install-hooks

# Verify all required local binaries
make check-env
```

---

## 4-Tier Local Verification Pyramid

Before committing or pushing, verify changes locally:

### Tier 1 — Static analysis, schema parity & Python unit tests (<10s, runs on `git commit`)

```bash
make lint
make test-unit
```

### Tier 2 — Nextflow DAG & `nf-test` stub verification (no containers required)

```bash
make test-stub
```

### Tier 3 — Containerized `nf-test` integration suite with miniature fixture DBs

```bash
make test
```

### Tier 4 — Local GitHub Actions workflow execution via `gh act`

```bash
make test-ci
# Or run the full pre-push gate:
make preflight
```

If working offline without Docker, set `SKIP_DOCKER_TESTS=1` when pushing with the pre-push hook enabled.

---

## Code style

- **Nextflow (`.nf`, `.config`)**: Validated with `nextflow lint .` (Nextflow >= 26.04.0).
- **Python (`bin/`, `tests/python/`, `tests/fixtures/`)**: Formatted and linted with `ruff` (line length 100,
  `Python >= 3.11`) and tested with `pytest`.
- **JSON, YAML, Markdown, and JS**: Formatted with [Prettier](https://prettier.io/).

---

## Branching and pull requests

- `main` is the stable release branch.
- Feature branches: `feature/<description>`
- Bug/hot fixes: `fix/<description>`
- Non-code tasks: `chore/<description>`
- Preparing a release: `release/<description>`
- Open a pull request against `main`; all CI checks must pass before merge.
