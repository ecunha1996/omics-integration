#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process createPathwaysMap {
    // Define input parameters
    input:
    val CONFIG_PATH

    output:
    val 'is_ok'

    script:
    """
    #!/usr/bin/env python
    import sys
    sys.path.append("$projectDir")
    from sampling import create_pathways_map
    create_pathways_map("$CONFIG_PATH")
    """
}

process sample{
    input:
    val ready
    val CONFIG_PATH

    output:
    val 'is_ok'

    script:
    """
    #!/usr/bin/env python
    import sys
    sys.path.append("$projectDir")
    from sampling import sample
    sample("$CONFIG_PATH")
    """
}

process samplingAnalysis{
    input:
    val ready
    val CONFIG_PATH

    output:
    val 'is_ok'

    script:
    """
    #!/usr/bin/env python
    import sys
    sys.path.append("$projectDir")
    from sampling import sampling_analysis
    sampling_analysis("$CONFIG_PATH")
    """
}

workflow {
         is_ok = createPathwaysMap(params.CONFIG_PATH)
         is_ok = sample(is_ok, params.CONFIG_PATH)
         samplingAnalysis(is_ok, params.CONFIG_PATH)
}