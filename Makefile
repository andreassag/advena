PIPELINE = main.nf
PROFILE ?= docker
ACT ?= gh act

.PHONY: help check-env install-hooks run resume lint test-unit test-stub test test-ci preflight clean

help:
	@echo "========================================================================"
	@echo "                  advena v2.0.0 Nextflow Pipeline                       "
	@echo "========================================================================"
	@echo "Available commands:"
	@echo "  make check-env     : Verify local dev, Docker, and gh-act binaries"
	@echo "  make run           : Run pipeline with default profile ($(PROFILE))"
	@echo "  make resume        : Resume previous run using Nextflow cache"
	@echo "  make lint          : Run Nextflow, Prettier, Ruff, schema & MkDocs checks"
	@echo "  make test-unit     : Run pytest unit tests on bin/ Python helpers"
	@echo "  make test-stub     : Run default & --run_all stub DAG + nf-test stubs"
	@echo "  make test          : Run full nf-test suite with $(PROFILE)"
	@echo "  make test-ci       : Run GitHub Actions PR workflows locally via gh act"
	@echo "  make preflight     : Run full 4-tier pre-push verification gate"
	@echo "  make clean         : Remove work, results, logs, and test caches"
	@echo "========================================================================"

check-env:
	@for cmd in nextflow nf-test python3 npx docker gh; do \
		command -v $$cmd >/dev/null 2>&1 || { echo "Missing required tool: $$cmd"; exit 1; }; \
	done
	@gh act --version >/dev/null 2>&1 || { echo "Missing gh-act extension. Install via: gh extension install https://github.com/nektos/gh-act"; exit 1; }
	@echo "All required local verification tools are available."

run:
	nextflow run $(PIPELINE) --input assets/samplesheet.csv --outdir results -profile $(PROFILE)

resume:
	nextflow run $(PIPELINE) --input assets/samplesheet.csv --outdir results -profile $(PROFILE) -resume

lint:
	nextflow lint .
	npx prettier --check .
	ruff check bin/ tests/python/ tests/fixtures/
	ruff format --check bin/ tests/python/ tests/fixtures/
	python3 tests/python/validate_schemas.py
	mkdocs build --strict

test-unit:
	pytest -v tests/python/

test-stub:
	nextflow run $(PIPELINE) -profile test -stub-run
	nextflow run $(PIPELINE) -profile test --run_all -stub-run
	nf-test test --tag stub --verbose

test:
	nf-test test --profile $(PROFILE) --verbose

test-ci:
	$(ACT) pull_request -W .github/workflows/linting.yml
	$(ACT) pull_request -W .github/workflows/schema-validate.yml
	$(ACT) pull_request -W .github/workflows/pipeline-test.yml
	$(ACT) pull_request -W .github/workflows/nf-test.yml -j nf-test
	$(ACT) pull_request -W .github/workflows/docs.yml -j build

preflight: check-env lint test-unit test-stub test test-ci

clean:
	@echo "Cleaning up Nextflow, nf-test, Python, and act artefacts..."
	rm -rf work/ results/ site/ .nextflow/ .nextflow.log .nextflow.log* .nf-test/ .act-artifacts/ .pytest_cache/ .ruff_cache/ .coverage
