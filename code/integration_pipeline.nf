#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Define the process to run your Python script
process getIntegrationResults {
    // Define input parameters
    input:
    val CONFIG_PATH
    script:
    """
    #!/usr/bin/env python
    import sys
    sys.path.append("$projectDir")
    from troppo_integration import integrate
    integrate("$CONFIG_PATH")
    """
}

process reconstructModels {
    // Define input parameters
    input:
    val CONFIG_PATH
    script:
    """
    #!/usr/bin/env python
    import sys
    sys.path.append("$projectDir")
    from troppo_integration import reconstruct_models
    reconstruct_models("$CONFIG_PATH")
    """
}


// Define the workflow
workflow {
    getIntegrationResults(params.CONFIG_PATH)
    reconstructModels(params.CONFIG_PATH)
}
