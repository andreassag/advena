#!/usr/bin/env nextflow

include { validateParameters; paramsHelp; paramsSummaryLog } from "plugin/nf-schema"

include { ADVENA } from "./workflows/advena"

//
// MAIN WORKFLOW
//
workflow {

    main:
    if (params.version) {
        log.info "${workflow.manifest.name} ${workflow.manifest.version}"
        exit 0
    }

    if (params.help) {
        def help_param = (params.help instanceof String && params.help != "true") ? params.help : ""
        log.info paramsHelp(
            command: "nextflow run ${workflow.manifest.name} --input samplesheet.csv --outdir results -profile docker",
            parameter: help_param
        )
        exit 0
    }

    validateParameters()
    log.info paramsSummaryLog()

    ADVENA()

    ADVENA.out.versions_blocks
        .collectFile(
            name: "software_versions.yml",
            storeDir: "${params.outdir}/pipeline_info",
            sort: true,
            newLine: false
        )

    onComplete:
    log.info(workflow.success ? "\nPipeline completed successfully!\n" : "Pipeline failed.")
}
