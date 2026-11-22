#!/usr/bin/env nextflow

include { validateParameters; paramsHelp; paramsSummaryLog } from "plugin/nf-schema"

include { ICESCREEN } from "./workflows/icescreen"

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

    ICESCREEN()

    channel
        .topic("versions")
        .distinct()
        .map { proc_name, tool, version ->
            def simple_name = proc_name.tokenize(":").last()
            [simple_name, "${tool}: \"${version}\""]
        }
        .groupTuple(by: 0)
        .map { simple_name, tool_entries ->
            def body = tool_entries.unique().sort().collect { entry -> "  ${entry}" }.join("\n")
            "${simple_name}:\n${body}\n"
        }
        .collectFile(
            name: "software_versions.yml",
            storeDir: "${params.outdir}/pipeline_info",
            sort: true,
            newLine: false
        )

    onComplete:
    log.info(workflow.success ? "\nPipeline completed successfully!\n" : "Pipeline failed.")
}
