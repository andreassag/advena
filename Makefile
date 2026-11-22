PIPELINE = main.nf
PROFILE ?= docker

.PHONY: help run resume test test-stub lint clean

help:
	@echo "========================================================================"
	@echo "                     advena Nextflow Pipeline                           "
	@echo "========================================================================"
	@echo "Available commands:"
	@echo "  make run        : Run the pipeline with default profile ($(PROFILE))"
	@echo "  make resume     : Resume the previous run using Nextflow cache"
	@echo "  make test-stub  : Run stub pipeline and stub nf-test suite"
	@echo "  make test       : Run full nf-test suite with Docker"
	@echo "  make lint       : Run Nextflow, Prettier, and schema checks"
	@echo "  make clean      : Remove work, results, logs, and test caches"
	@echo "========================================================================"

run:
	nextflow run $(PIPELINE) --input assets/samplesheet.csv --outdir results -profile $(PROFILE)

resume:
	nextflow run $(PIPELINE) --input assets/samplesheet.csv --outdir results -profile $(PROFILE) -resume

test-stub:
	nextflow run $(PIPELINE) -profile test -stub-run
	nf-test test --tag stub --verbose

test:
	nf-test test --profile docker --verbose

lint:
	nextflow lint .
	npx prettier --check .

clean:
	@echo "Cleaning up Nextflow and nf-test artefacts..."
	rm -rf work/ results/ site/ .nextflow/ .nextflow.log .nextflow.log* .nf-test/
