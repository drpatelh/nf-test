import org.yaml.snakeyaml.Yaml

// Adapted from nf-core/rnaseq source code for lib/WorkflowMain.groovy. 
// See CITATIONS.md for link to nf-core

class Helper {
    //
    // Print parameter summary log to screen
    //
    public static String paramsSummaryLog(workflow, params, log) {
        def summary_log = ''
        summary_log += SchemaValidation.paramsSummaryLog(workflow, params)
        summary_log += Logging.dashedLine()
        return summary_log
    }

    public static Map getqcAndArchRParams(paramPath){
        String configText = new File(paramPath).text
        def config = new Yaml().load(configText) 
        return config
    }
    //
    // Citation string for pipeline
    //
    public static String citation(workflow) {
        return "If you use ${workflow.manifest.name} for your analysis please cite:\n\n" +
            "Dummy info here"
    }

    //
    // Print help to screen if required
    //
    public static String help(workflow, params, log) {
        def command = "nextflow run ${workflow.manifest.name} -profile docker -params-file <params.yml> --outDir <OUTDIR> "
        def help_string = ''
        help_string += SchemaValidation.paramsHelp(workflow, params, command)
        help_string += Logging.dashedLine()
        return help_string
    }

    //
    // Validate parameters and print summary to screen
    //
    public static void initialise(workflow, params, log) {
        // Write logo 
        print Logging.logo(workflow)
        if (params.help) {
            print help(workflow, params, log)
            System.exit(0)
        }

        // Validate workflow parameters via the JSON schema
        SchemaValidation.validateAdditionalParameters(workflow, params, log)
        SchemaValidation.validateWorkflowParameters(workflow, params, log)

        // Print parameter summary log to screen
        log.info paramsSummaryLog(workflow, params, log)

        // Check that a -profile or Nextflow config has been provided to run the pipeline
        Logging.checkConfigProvided(workflow, log)

        // Check AWS batch settings
        Logging.awsBatch(workflow, params)
    }
}

