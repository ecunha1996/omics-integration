import json
import math
import os
from os.path import abspath, dirname, join
import pandas as pd
from gsmmutils import DATA_PATH
from gsmmutils.io import read_csv
from gsmmutils.omics.omics_integration import OmicsIntegration
from gsmmutils.omics.omics_processing import thresholding_filter
from gsmmutils.omics.troppo_integration import troppo_omics_integration
from gsmmutils.omics.model_handle import load_model, sbml_model_reconstruction
import matplotlib.pyplot as plt
import argparse
from numpy import linspace

def reconstruction_pipeline():
    """
    This function is used to run the reconstruction pipeline.

    In the end this function generates a SBML file of the tissue-specific reconstructed_models that resulted from the
    omics data integration with Troppo.

    All the parameters required to run this pipeline can be defined in the ***pipeline_paths.py***,
    ***config_variables.py***, and ***medium_variables.py files***.
    """
    print('-------------------------------------------------------------------------------------------------------')
    print('--------------------------------------- Loading Template Model. ---------------------------------------')
    print('-------------------------------------------------------------------------------------------------------')

    template_model = load_model(model_path=MODEL_PATH, consistent_model_path=CONSISTENT_MODEL_PATH, MEDIUM_CONDITIONS=MEDIUM_CONDITIONS)
    print('-------------------------------------------------------------------------------------------------------')
    print('-------------------------------------- Processing Omics Dataset. --------------------------------------')
    print('-------------------------------------------------------------------------------------------------------')

    omics = OmicsIntegration(OMICS_DATA_PATH, samples_names={
    }, model=template_model)

    omics.getmm = read_csv(OMICS_DATA_PATH, index_name='GeneID', index_col=0, comment='#', sep='\t')
    omics.sum_tech_reps()
    omics_data = omics.counts.applymap(lambda x: math.log2(x + 1))
    print(omics_data.describe())
    fig = omics_data.plot.density()
    fig.set_xticks([x / 2 for x in range(0, 40)])
    plt.savefig(rf"{DATA_PATH}/ngaditana/PRJNA589063/counts_density_light.png")
    omics_data = omics_data.T
    print('Omics dataset Loaded.')

    if params['THRESHOLDING_STRATEGY'] != 'default':
        omics_data = thresholding_filter(omics_dataframe=omics_data,
                                         thresholding_strategy=params['THRESHOLDING_STRATEGY'],
                                         global_threshold_upper=params['GLOBAL_THRESHOLD_UPPER'],
                                         global_threshold_lower=params['GLOBAL_THRESHOLD_LOWER'],
                                         local_threshold=params['LOCAL_THRESHOLD'])

        print(f'{params["THRESHOLDING_STRATEGY"]} threshold filter applied.')

    print('-------------------------------------------------------------------------------------------------------')
    print('------------------------------- Starting Omics Integration with Troppo. -------------------------------')
    print('-------------------------------------------------------------------------------------------------------')

    integration_result = {}

    for algorithm in params["ALGORITHMS"]:
        troppo_result_dict = {}
        print(f'Omics integration with {algorithm} started.')
        print('-------------------------------------------------------------------------------------------------------')

        if algorithm == 'fastcore':
            thred_number = params["THREAD_NUMBER_FASTCORE"]

        elif algorithm == 'tinit':
            thred_number = params["THREAD_NUMBER_TINIT"]

        else:
            return 'Algorithm not supported by this pipeline.'

        for threshold in INTEGRATION_THRESHOLDS:
            troppo_result = troppo_omics_integration(model=template_model, algorithm=algorithm, threshold=threshold,
                                                     thread_number=thred_number, omics_dataset=omics_data,
                                                     params=params)

            for sample in list(troppo_result.keys()):
                th = str(round(threshold, 2)).replace('.', '_')
                troppo_result_dict[f'{params["MODEL"]}_{sample}_{algorithm}_t{th}'] = troppo_result[sample]
                integration_result[f'{params["MODEL"]}_{sample}_{algorithm}_t{th}'] = troppo_result[sample]

            print('----------------------------------------------------------'
                  '---------------------------------------------')

        if params["THRESHOLDING_STRATEGY"] == 'default':
            result_path = os.path.join(TROPPO_RESULTS_PATH,
                                       f'{params["MODEL"]}_{params["DATASET"]}_{algorithm}_{params["THRESHOLDING_STRATEGY"]}.csv')
        else:
            file_path = f'{params["MODEL"]}_{params["DATASET"]}_{algorithm}_{params["THRESHOLDING_STRATEGY"]}_' \
                        f'{params["GLOBAL_THRESHOLD_UPPER"]}_{params["GLOBAL_THRESHOLD_LOWER"]}_{params["LOCAL_THRESHOLD"]}.csv'
            result_path = os.path.join(TROPPO_RESULTS_PATH, file_path)

        troppo_result_dataframe = pd.DataFrame.from_dict(troppo_result_dict, orient='index')
        troppo_result_dataframe.to_csv(result_path)

    if params["EVALUATE_TASKS"]:
        print('------------------------ Starting task evaluation for the integration results. ------------------------')
        print('-------------------------------------------------------------------------------------------------------')

        # task_eval(model_template=template_model, integration_results=integration_result)

        print('-------------------------------------------------------------------------------------------------------')

    if params["RECONSTRUCT_MODELS"]:
        print('----------------------- Starting reconstruction of the context-specific models. -----------------------')
        print('-------------------------------------------------------------------------------------------------------')

        for sample_name in list(integration_result.keys()):
            print(f'Context-specific model reconstruction for {sample_name} started.')
            print('---------------------------------------------------------------'
                  '----------------------------------------')
            sbml_model_reconstruction(original_model=template_model, sample=sample_name,
                                      integration_result_dict=integration_result,
                                      MEDIUM_CONDITIONS=MEDIUM_CONDITIONS)

            print('---------------------------------------------------------------'
                  '----------------------------------------')

    print('------------------------------------------ Pipeline Finished ------------------------------------------')
    print('-------------------------------------------------------------------------------------------------------')


def read_args():
    parser = argparse.ArgumentParser(description='Process some integers.')

    # Add the arguments
    parser.add_argument('--CONFIG', type=str, required=False, help='Configuration file')
    parser.add_argument('--DATASET', type=str, required=False, help='The dataset to use')
    parser.add_argument('--MODEL', type=str, required=False, help='The model to use')
    parser.add_argument('--RAW_DATA', type=str, required=False, help='The raw data to use')
    parser.add_argument('--OBJECTIVE', type=str, required=False, help='The objective to use')
    parser.add_argument('--MEDIUM_NAME', type=str, required=False, help='The medium name to use')
    parser.add_argument('--NOMENCLATURE', type=str, required=False, help='The nomenclature to use')
    parser.add_argument('--OMICS_TYPE', type=str, required=False, help='The omics type to use')
    parser.add_argument('--THRESHOLDING_STRATEGY', type=str, required=False, help='The thresholding strategy to use')
    parser.add_argument('--GLOBAL_THRESHOLD_UPPER', type=int, required=False, help='The global threshold upper to use')
    parser.add_argument('--GLOBAL_THRESHOLD_LOWER', type=int, required=False, help='The global threshold lower to use')
    parser.add_argument('--LOCAL_THRESHOLD', type=int, required=False, help='The local threshold to use')
    parser.add_argument('--ALGORITHMS', type=str, nargs='+', required=False, help='The algorithms to use')
    parser.add_argument('--THREAD_NUMBER_FASTCORE', type=int, required=False, help='The thread number for fastcore')
    parser.add_argument('--THREAD_NUMBER_TINIT', type=int, required=False, help='The thread number for tinit')
    parser.add_argument('--EVALUATE_TASKS', type=str, required=False, help='Whether to evaluate tasks')
    parser.add_argument('--GAP_FILLING', type=str, required=False, help='Whether to perform gap filling')
    parser.add_argument('--RECONSTRUCT_MODELS', type=str, required=False, help='Whether to reconstruct models')

    args = parser.parse_args()

    return vars(args)

if __name__ == '__main__':
    input_params = read_args()
    DATA_PATH = "/home/ecunha/omics-integration/data"
    os.chdir(DATA_PATH)
    CONFIG_PATH = input_params.get("CONFIG") or abspath(join(dirname(__file__), '../config'))
    params = json.load(open(rf"{CONFIG_PATH}/troppo/PRJNA589063.json", "r"))
    UPTAKE_DRAINS = {
        'f2 medium': {"EX_C00014__dra", "EX_C00059__dra", "EX_C01330__dra", "EX_C00011__dra", "EX_C14818__dra",
                      "EX_C00080__dra", "EX_C00001__dra", "EX_C00305__dra", "EX_C01327__dra", "EX_C00244__dra",
                      "EX_C00009__dra",
                      "EX_C00007__dra", "EX_C00205__dra", "EX_C00378__dra", "EX_C00120__dra", "EX_C02823__dra",
                      'DM_C00369__chlo'
                      }}
    params['uptake_drains'] = UPTAKE_DRAINS
    MODEL_PATH = input_params.get("MODEL_PATH") or params.get("MODEL_PATH")
    CONSISTENT_MODEL_PATH = input_params.get("CONSISTENT_MODEL_PATH") or params.get("CONSISTENT_MODEL_PATH")
    OMICS_DATA_PATH = input_params.get("RAW_DATA") or params.get("RAW_DATA")
    TROPPO_RESULTS_PATH = f"{DATA_PATH}/integration/{params['DATASET']}"
    MODEL_RESULTS_PATH = f"{DATA_PATH}/integration/{params['DATASET']}"
    MODEL_TASKS_PATH = f"{DATA_PATH}/integration"
    MODEL_TASKS_RESULTS_PATH = f"{DATA_PATH}/integration/{params['DATASET']}"
    media = pd.read_excel(rf"{DATA_PATH}/ngaditana/media.xlsx", index_col=0, sheet_name=None, engine='openpyxl')[
        'media_with_starch'].to_dict(orient='index')
    MEDIUM_CONDITIONS = {'f2 medium': {key: (value['LB'], value['UB']) for key, value in media.items()}}

    MEDIUM_METABOLITES = {'f2 medium': [e.split("__")[0].replace("EX_", "") for e in MEDIUM_CONDITIONS['f2 medium']]}

    # THRESHOLDS = linspace(0, 9, 10)  # List of ints or floats.
    THRESHOLDS = [2]
    INTEGRATION_THRESHOLDS = list(THRESHOLDS)
    reconstruction_pipeline()
    print("")
