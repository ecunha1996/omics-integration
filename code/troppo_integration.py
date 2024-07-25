import json
import math
import os
from cgitb import reset
from os.path import abspath, dirname, join

import numpy as np
import pandas as pd
from gsmmutils import DATA_PATH
from gsmmutils.io import read_csv
from gsmmutils.omics.omics_integration import OmicsIntegration
from gsmmutils.omics.omics_processing import thresholding_filter
from gsmmutils.omics.troppo_integration import troppo_omics_integration
from gsmmutils.omics.model_handle import load_model, sbml_model_reconstruction
import matplotlib.pyplot as plt
import argparse
from pca_analysis import rank_models_difference
from numpy import linspace
from scipy.interpolate import interp1d
from troppo.omics import GeneLevelThresholding


def integrate(CONFIG_PATH):
    params = json.load(open(rf"{CONFIG_PATH}", "r"))
    media_path = params.get("MEDIA_PATH")
    media_sheet = params.get("MEDIA_SHEET")
    model_path = params.get("MODEL_PATH")
    consistent_model_path  = params.get("CONSISTENT_MODEL_PATH")
    raw_data_path =  params.get("RAW_DATA")
    getmm_path = params.get("GETMM")
    media = pd.read_excel(media_path, index_col=0, sheet_name=None, engine='openpyxl')[media_sheet].to_dict(orient='index')
    medium_conditions = {'f2 medium': {key: (value['LB'], value['UB']) for key, value in media.items()}}
    template_model = load_model(model_path=model_path, consistent_model_path=consistent_model_path, medium_conditions=medium_conditions)
    omics = OmicsIntegration(raw_data_path, samples_names={
    }, model=template_model)
    if params.get("SAMPLES_TO_DROP"):
        for sample in params.get("SAMPLES_TO_DROP"):
            omics.drop_sample(sample)
    omics.getmm = pd.read_csv(getmm_path, index_col=0, sep="\t")
    omics.sum_tech_reps()
    omics_data = omics.getmm.applymap(lambda x: math.log2(x + 1))
    omics_data.index = [e.replace("-", "_").replace("_1", "") for e in omics_data.index]
    omics_data = omics_data.loc[omics_data.index.isin([gene.id for gene in template_model.genes])]
    omics_data = omics_data.T
    # thresholds = {0.3: 0.3}
    print('Omics dataset Loaded.')
    threshold = GeneLevelThresholding(omics_dataframe=omics_data,
                                      thresholding_strat=params['THRESHOLDING_STRATEGY'],
                                      global_threshold_lower=params['GLOBAL_THRESHOLD_LOWER'],
                                      global_threshold_upper=params['GLOBAL_THRESHOLD_UPPER'],
                                      local_threshold=params['LOCAL_THRESHOLD'])
    omics_data = threshold.apply_thresholding_filter()
    thresholds = {round(quantile, 3): np.quantile(omics_data, quantile).round(3) for quantile in linspace(0.10, 0.35, 3)} #51
    print(f"Thresholds: {thresholds}")
    for algorithm in params["ALGORITHMS"]:
        troppo_result_dict = {}
        thread_number = params["THREAD_NUMBER"]
        for quantile, threshold in thresholds.items():
            troppo_result = troppo_omics_integration(model=template_model, algorithm=algorithm, threshold=threshold,
                                                     thread_number=thread_number, omics_dataset=omics_data,
                                                     params=params)

            for sample in list(troppo_result.keys()):
                th = str(round(quantile, 3))
                troppo_result_dict[f'{sample.split("_")[0]}_{algorithm}_{th}'] = troppo_result[sample]

        file_path = f'{algorithm}_{params["THRESHOLDING_STRATEGY"].replace(" ", "_")}_' \
                    f'{params["GLOBAL_THRESHOLD_UPPER"]}_{params["GLOBAL_THRESHOLD_LOWER"]}_{params["LOCAL_THRESHOLD"]}.csv'
        result_path = os.path.join(params.get("TROPPO_RESULTS_PATH"), 'integration', file_path)
        troppo_result_dataframe = pd.DataFrame.from_dict(troppo_result_dict, orient='index')
        troppo_result_dataframe.to_csv(result_path)


def get_best_thresholds(df):
    dist = rank_models_difference(df)
    sorted_thresholds = sorted(dist, key=lambda k: dist[k], reverse=True)
    return sorted_thresholds[0]

def reconstruct_models(CONFIG_PATH):
    params = json.load(open(rf"{CONFIG_PATH}", "r"))
    media_path = params.get("MEDIA_PATH")
    media_sheet = params.get("MEDIA_SHEET")
    model_path = params.get("MODEL_PATH")
    consistent_model_path = params.get("CONSISTENT_MODEL_PATH")
    media = pd.read_excel(media_path, index_col=0, sheet_name=None, engine='openpyxl')[media_sheet].to_dict(orient='index')
    medium_conditions = {'f2 medium': {key: (value['LB'], value['UB']) for key, value in media.items()}}
    params = json.load(open(rf"{CONFIG_PATH}", "r"))
    template_model = load_model(model_path=model_path, consistent_model_path=consistent_model_path, medium_conditions=medium_conditions)
    for algorithm in params.get("ALGORITHMS"):
        file_path = f'{algorithm}_{params["THRESHOLDING_STRATEGY"].replace(" ", "_")}_' \
                    f'{params["GLOBAL_THRESHOLD_UPPER"]}_{params["GLOBAL_THRESHOLD_LOWER"]}_{params["LOCAL_THRESHOLD"]}.csv'
        df = pd.read_csv(os.path.join(params.get("TROPPO_RESULTS_PATH"), 'integration',file_path), index_col=0)
        df = df.astype(int)
        best_threshold =get_best_thresholds(df)
        integration_result = {}
        df = df.filter(regex=f".*{best_threshold}$", axis=0)
        for row in df.iterrows():
            integration_result[row[0]] = row[1]
        for sample_name in list(integration_result.keys()):
            sbml_model_reconstruction(original_model=template_model, sample=str(sample_name),
                                      integration_result_dict=integration_result, params=params,
                                      medium_conditions=medium_conditions,
                                      model_results_path=os.path.join(params.get("TROPPO_RESULTS_PATH"), 'integration/models'))


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

    template_model = load_model(model_path=MODEL_PATH, consistent_model_path=consistent_model_path, medium_conditions=medium_conditions)

    print('-------------------------------------------------------------------------------------------------------')
    print('-------------------------------------- Processing Omics Dataset. --------------------------------------')
    print('-------------------------------------------------------------------------------------------------------')

    omics = OmicsIntegration(OMICS_DATA_PATH, samples_names={
    }, model=template_model)
    omics.drop_sample("C3_4")
    omics.getmm = pd.read_csv(GETMM_PATH, index_col=0, sep="\t")
    omics.sum_tech_reps()
    omics_data = omics.getmm.applymap(lambda x: math.log2(x + 1))
    print(omics_data.describe())
    omics_data.index = [e.replace("-", "_").replace("_1", "") for e in omics_data.index]
    omics_data = omics_data.loc[omics_data.index.isin([gene.id for gene in template_model.genes])]
    # fig = omics_data.plot.density()
    # plt.show()
    omics_data = omics_data.T
    # for condition in omics_data.index:
    thresholds = {0.3: 0.3}

    # thresholds = {quantile.round(2): quantile.round(2) for quantile in linspace(1, 5, 401)}
    #
    print('Omics dataset Loaded.')

    if params['THRESHOLDING_STRATEGY'] != 'default':
        threshold = GeneLevelThresholding(omics_dataframe=omics_data,
                                          thresholding_strat=params['THRESHOLDING_STRATEGY'],
                                          global_threshold_lower=params['GLOBAL_THRESHOLD_LOWER'],
                                          global_threshold_upper=params['GLOBAL_THRESHOLD_UPPER'],
                                          local_threshold=params['LOCAL_THRESHOLD'])
        omics_data = threshold.apply_thresholding_filter()

        print(f'{params["THRESHOLDING_STRATEGY"]} threshold filter applied.')
    print('-------------------------------------------------------------------------------------------------------')
    print('------------------------------- Starting Omics Integration with Troppo. -------------------------------')
    print('-------------------------------------------------------------------------------------------------------')
    # thresholds = {round(quantile, 3): np.quantile(omics_data, quantile).round(3) for quantile in linspace(0.10, 0.35, 51)}
    print(thresholds)
    integration_result = {}

    for algorithm in params["ALGORITHMS"]:
        troppo_result_dict = {}
        print(f'Omics integration with {algorithm} started.')
        print('-------------------------------------------------------------------------------------------------------')

        thred_number = params["THREAD_NUMBER"]
        # for condition, cond_thresholds in thresholds.items():
        #     sample_name = f"{condition}_{params['THRESHOLDING_STRATEGY']}_{params['GLOBAL_THRESHOLD_UPPER']}_{params['GLOBAL_THRESHOLD_LOWER']}_{params['LOCAL_THRESHOLD']}"
        for quantile, threshold in thresholds.items():
            troppo_result = troppo_omics_integration(model=template_model, algorithm=algorithm, threshold=threshold,
                                                     thread_number=thred_number, omics_dataset=omics_data,
                                                     params=params)

            for sample in list(troppo_result.keys()):
                th = str(round(quantile, 3))
                troppo_result_dict[f'{sample.split("_")[0]}_{algorithm}_{th}'] = troppo_result[sample]
                integration_result[f'{sample.split("_")[0]}_{algorithm}_{th}'] = troppo_result[sample]

            print('----------------------------------------------------------'
                  '---------------------------------------------')

        if params["THRESHOLDING_STRATEGY"] == 'default':
            result_path = os.path.join(TROPPO_RESULTS_PATH,
                                       f'{algorithm}_{params["THRESHOLDING_STRATEGY"]}.csv')
        else:
            file_path = f'{algorithm}_{params["THRESHOLDING_STRATEGY"].replace(" ", "_")}_' \
                        f'{params["GLOBAL_THRESHOLD_UPPER"]}_{params["GLOBAL_THRESHOLD_LOWER"]}_{params["LOCAL_THRESHOLD"]}.csv'
            result_path = os.path.join(TROPPO_RESULTS_PATH, file_path)

        troppo_result_dataframe = pd.DataFrame.from_dict(troppo_result_dict, orient='index')
        troppo_result_dataframe.to_csv(result_path)

        # df =  pd.read_csv(os.path.join(TROPPO_RESULTS_PATH, rf'{algorithm}_{params["THRESHOLDING_STRATEGY"].replace(" ", "_")}_'f'{params["GLOBAL_THRESHOLD_UPPER"]}_{params["GLOBAL_THRESHOLD_LOWER"]}_{params["LOCAL_THRESHOLD"]}.csv'), index_col=0)
        # df = df.filter(regex=f".*1.2$", axis=0)
        # for row in df.iterrows():
        #     integration_result[row[0]] = row[1]

    if params["RECONSTRUCT_MODELS"]:
        print('----------------------- Starting reconstruction of the context-specific models. -----------------------')
        print('-------------------------------------------------------------------------------------------------------')

        for sample_name in list(integration_result.keys()):
            print(f'Context-specific model reconstruction for {sample_name} started.')
            print('---------------------------------------------------------------'
                  '----------------------------------------')
            sbml_model_reconstruction(original_model=template_model, sample=str(sample_name),
                                      integration_result_dict=integration_result, params=params,
                                      medium_conditions=medium_conditions,
                                      model_results_path= MODEL_RESULTS_PATH)

            print('---------------------------------------------------------------'
                  '----------------------------------------')

    print('------------------------------------------ Pipeline Finished ------------------------------------------')
    print('-------------------------------------------------------------------------------------------------------')


def read_args():
    parser = argparse.ArgumentParser(description='Process some integers.')

    # Arguments to add
    args_list = [
        ('--CONFIG', str, 'Configuration file'),
        ('--DATASET', str, 'The dataset to use'),
        ('--MODEL', str, 'The model to use'),
        ('--RAW_DATA', str, 'The raw data to use'),
        ('--OBJECTIVE', str, 'The objective to use'),
        ('--MEDIUM_NAME', str, 'The medium name to use'),
        ('--NOMENCLATURE', str, 'The nomenclature to use'),
        ('--OMICS_TYPE', str, 'The omics type to use'),
        ('--THRESHOLDING_STRATEGY', str, 'The thresholding strategy to use'),
        ('--GLOBAL_THRESHOLD_UPPER', int, 'The global threshold upper to use'),
        ('--GLOBAL_THRESHOLD_LOWER', int, 'The global threshold lower to use'),
        ('--LOCAL_THRESHOLD', int, 'The local threshold to use'),
        ('--ALGORITHMS', str, 'The algorithms to use', '+'),
        ('--THREAD_NUMBER_FASTCORE', int, 'The thread number for fastcore'),
        ('--THREAD_NUMBER_TINIT', int, 'The thread number for tinit'),
        ('--GAP_FILLING', str, 'Whether to perform gap filling'),
        ('--RECONSTRUCT_MODELS', str, 'Whether to reconstruct models')
    ]

    for arg, arg_type, help_text, *nargs in args_list:
        parser.add_argument(arg, type=arg_type, required=False, help=help_text, nargs=nargs[0] if nargs else None)

    return vars(parser.parse_args())



# if __name__ == '__main__':
#     input_params = read_args()
#     DATA_PATH = "/home/ecunha/omics-integration/data"
#     RESULTS_PATH = "/home/ecunha/omics-integration/results"
#     os.chdir(DATA_PATH)
#     CONFIG_PATH = input_params.get("CONFIG") or abspath(join(dirname(__file__), '../config'))
#     params = json.load(open(rf"{CONFIG_PATH}", "r"))
#     UPTAKE_DRAINS = {
#         'f2 medium': {"EX_C00014__dra", "EX_C00059__dra", "EX_C01330__dra", "EX_C00011__dra", "EX_C14818__dra",
#                       "EX_C00080__dra", "EX_C00001__dra", "EX_C00305__dra", "EX_C01327__dra", "EX_C00244__dra",
#                       "EX_C00009__dra",
#                       "EX_C00007__dra", "EX_C00205__dra", "EX_C00378__dra", "EX_C00120__dra", "EX_C02823__dra",
#                       'DM_C00965__chlo'
#                       }}
#     params['uptake_drains'] = UPTAKE_DRAINS
#     MODEL_PATH = input_params.get("MODEL_PATH") or params.get("MODEL_PATH")
#     consistent_model_path = input_params.get("CONSISTENT_MODEL_PATH") or params.get("CONSISTENT_MODEL_PATH")
#     OMICS_DATA_PATH = input_params.get("RAW_DATA") or params.get("RAW_DATA")
#     GETMM_PATH = input_params.get("GETMM") or params.get("GETMM")
#     TROPPO_RESULTS_PATH = f"{RESULTS_PATH}/ngaditana/{params['DATASET']}/integration"
#     MODEL_RESULTS_PATH = join(TROPPO_RESULTS_PATH, "models")
#
#     MEDIUM_METABOLITES = {'f2 medium': [e.split("__")[0].replace("EX_", "") for e in medium_conditions['f2 medium']]}
#     reconstruction_pipeline()
#     print("")