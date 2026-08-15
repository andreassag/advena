# Contributing

Contributions are welcome. This document describes how to set up a development environment and the conventions the
project follows.

---

## Development setup

```bash
git clone https://github.com/andreassag/advena.git
cd advena

# Install Node.js dependencies (Prettier)
npm install

# Install and activate the pre-commit hooks
pip install pre-commit
pre-commit install
```

---

## Code style

### Nextflow (`.nf`, `.config`)

Nextflow files are formatted with [Prettier](https://prettier.io/) via `prettier-plugin-groovy`.

---

## Pre-commit hooks

The following hooks run automatically on `git commit`:

| Hook                  | Purpose                                       |
| --------------------- | --------------------------------------------- |
| `prettier`            | Formatting for Nextflow, JSON, YAML, Markdown |
| `trailing-whitespace` | Strip trailing whitespace                     |
| `end-of-file-fixer`   | Ensure files end with a newline               |
| `check-yaml`          | Validate YAML syntax                          |
| `check-json`          | Validate JSON syntax                          |
| `mixed-line-ending`   | Enforce LF line endings                       |

Run all hooks manually:

```bash
pre-commit run --all-files
```

---

## Testing

### Stub run (fast, no software required)

```bash
nextflow run . -profile conda,test -stub-run
```

All 34 processes should complete successfully. This validates the channel topology without executing any real tools.

### Full pipeline test

```bash
nextflow run . \
    --input assets/samplesheet.csv \
    --outdir results \
    -profile test,docker
```

Requires Docker and an accessible ICEscreen installation (`ICESCREEN_ROOT`).

### nf-test

```bash
# requires Nextflow and Docker
nf-test test --profile docker --verbose
```

Tests live in `tests/`.

---

## Branching and pull requests

- `main` is the stable release branch.
- Feature branches: `feature/<description>`
- Bug/hot fixes: `fix/<description>`
- Non-code tasks: `chore/<description>`
- Preparing a release: `release/<description>`
- Open a pull request against `main`; CI must pass before merge.

---

## Releases

Releases follow [Semantic Versioning](https://semver.org/). To cut a release:

1. Update the version in `nextflow.config` (`manifest.version`).
2. Commit with message `release: bump version to X.Y.Z`.
3. Push a `X.Y.Z` tag — the `release.yml` workflow creates the GitHub Release automatically.

```bash
git tag 1.1.0
git push origin 1.1.0
```
