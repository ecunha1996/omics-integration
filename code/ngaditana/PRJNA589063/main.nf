#!/usr/bin/env nextflow

params.model = "../data/ngaditana/models/model_ng.xml"
params.specific_models_dir = "../results/ngaditana/PRJNA589063/integration/models"
params.results_dir = "../results/ngaditana/PRJNA589063/dfa/"

process generate_pathway_map {
    input:
    path modelFile

    output:
    path "pathways_map.csv"

    script:
    """
    python -c "
import os
from gsmmutils.model.COBRAmodel import MyModel

model = MyModel('${modelFile}', 'e_Biomass__cytop')
model.get_pathway_reactions_map()
results_dataframe = pd.DataFrame.from_dict(data=model.pathway_reactions_map, orient='index').T
results_dataframe.to_csv('pathways_map.csv', index=False)
"
    """
}

process achr_sampling {
    input:
    path modelFile

    output:
    path "${modelFile.simpleName}_ACHR_samples.csv" into sampled_files

    script:
    """
    python -c "
from sampling import achr_sample
import pandas as pd

filename = '${modelFile}'
biomass_reaction = 'e_Biomass__cytop'
samples = achr_sample(filename, biomass_reaction)
sample_df = pd.DataFrame(samples[f"{filename.split('/')[-1].split('.xml')[0]}"])
sample_df.to_csv(f"{filename.split('/')[-1].split('.xml')[0]}_ACHR_samples.csv", index=False, sep='\t')
"
    """
}

process kstest_analysis {
    input:
    path sampledFile1, sampledFile2 from sampled_files.collect().combinations().flatten().unique().filter { it[0].name.split('_')[0] == it[1].name.split('_')[0] && it[0].name.split('_')[1] == it[1].name.split('_')[1] }

    output:
    path "*_all_results.tsv" into kstest_results

    script:
    """
    python -c "
from sampling import kstest
import pandas as pd

samples_control = pd.read_csv('${sampledFile1}', sep='\t')
samples_condition = pd.read_csv('${sampledFile2}', sep='\t')
dataset_name = '${sampledFile1.simpleName}_${sampledFile2.simpleName}'
results_dir = '${params.results_dir}'
kstest(samples_control, samples_condition, dataset_name, results_dir)
"
    """
}

process pathway_enrichment {
    input:
    path kstestResult from kstest_results
    path pathwaysMap from pathways_map

    output:
    path "*_results.csv" into enrichment_results

    script:
    """
    python -c "
from sampling import pathway_enrichment
import pandas as pd

kstest_result = pd.read_csv('${kstestResult}', sep='\t')
rxnlist = kstest_result[kstest_result['Reject'] == True]['Reaction'].to_list()
dataset_name = '${kstestResult.simpleName}'
pathways_data = pd.read_csv('${pathwaysMap}')
results_dir = '${params.results_dir}'
pathway_enrichment(rxnlist, dataset_name, pathways_data, results_dir)
"
    """
}

workflow {
    generate_pathway_map(params.model)

    achr_sampling("${params.specific_models_dir}/*.xml")

    kstest_analysis(sampled_files.collect().combinations().flatten().unique().filter { it[0].name.split('_')[0] == it[1].name.split('_')[0] && it[0].name.split('_')[1] == it[1].name.split('_')[1] })

    pathway_enrichment(kstest_results, pathways_map)
}
