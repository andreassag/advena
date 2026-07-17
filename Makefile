# Variables - Customize these for your pipeline
PIPELINE = main.nf
PROFILE = conda

# Default target when you just run 'make'
.PHONY: help
help:
	@echo "========================================================================"
	@echo "                     Nextflow Pipeline Developer Tool                   "
	@echo "========================================================================"
	@echo "Available commands:"
	@echo "  make run      : Run the pipeline with default profile/configs"
	@echo "  make resume   : Resume the previous run using Nextflow cache"
	@echo "  make clean    : Remove work directory, execution logs, and cache"
	@echo "========================================================================"

# Run the pipeline from scratch
.PHONY: run
run:
	nextflow run $(PIPELINE) --input assets/samplesheet.csv --outdir results -profile $(PROFILE)

# Resume a stalled or tweaked pipeline run using cached tasks
.PHONY: resume
resume:
	nextflow run $(PIPELINE) --input assets/samplesheet.csv --outdir results -profile $(PROFILE) -resume

# Standard cleanup of temporary/working files
.PHONY: clean
clean:
	@echo "Cleaning up local Nextflow work and log tracks..."
	rm -rf work/ .nextflow/ .nextflow.log .nextflow.log*
