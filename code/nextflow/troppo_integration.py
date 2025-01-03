import json
import math
import os
import numpy as np
import pandas as pd
from cobra.flux_analysis import pfba
from gsmmutils.omics.model_handle import load_model, sbml_model_reconstruction
from gsmmutils.omics.omics_integration import OmicsIntegration
from gsmmutils.omics.troppo_integration import troppo_omics_integration
from numpy import linspace
from troppo.omics import GeneLevelThresholding
from pca_analysis import rank_models_difference


def preprocessing(params: dict):
    media = pd.read_excel(params.get("MEDIA_PATH"), index_col=0, sheet_name=None, engine='openpyxl')[params.get("MEDIA_SHEET")].to_dict(orient='index')
    medium_conditions = {'f2 medium': {key: (value['LB'], value['UB']) for key, value in media.items()}}
    # cobra.core.Configuration().solver = "cplex"
    template_model = load_model(model_path=params.get("MODEL_PATH"), consistent_model_path=params.get("CONSISTENT_MODEL_PATH"), medium_conditions=medium_conditions)
    omics = OmicsIntegration(params.get("RAW_DATA"), samples_names={}, model=template_model)
    if params.get("SAMPLES_TO_DROP"):
        for sample in params.get("SAMPLES_TO_DROP"):
            omics.drop_sample(sample)
    omics.getmm = pd.read_csv(params.get("GETMM"), index_col=0, sep="\t")
    omics.sum_tech_reps()
    omics_data = omics.getmm.applymap(lambda x: math.log2(x + 1))
    omics_data = omics_data.loc[omics_data.index.isin([gene.id for gene in template_model.genes])]
    omics_data = omics_data.T
    return omics_data, template_model

def integrate(config_path: str):
    params = json.load(open(rf"{config_path}", "r"))
    omics_data, template_model = preprocessing(params)
    print('Omics dataset Loaded.')
    threshold = GeneLevelThresholding(omics_dataframe=omics_data,
                                      thresholding_strat=params['THRESHOLDING_STRATEGY'],
                                      global_threshold_lower=params['GLOBAL_THRESHOLD_LOWER'],
                                      global_threshold_upper=params['GLOBAL_THRESHOLD_UPPER'],
                                      local_threshold=params['LOCAL_THRESHOLD'])
    omics_data = threshold.apply_thresholding_filter()
    # thresholds = {round(quantile, 3): np.quantile(omics_data, quantile).round(3) for quantile in linspace(0.05, 0.95, 19)}
    #thresholds = {round(quantile, 3): 00965np.quantile(omics_data, quantile).round(3) for quantile in linspace(0.4, 0.45, 1)}
    thresholds = {0.4: 0.40}
    print(f"Thresholds: {thresholds}")
    for algorithm in params["ALGORITHMS"]:
        troppo_result_dict = {}
        thread_number = params["THREAD_NUMBER"]
        for quantile, threshold in thresholds.items():
            omics_data = omics_data.loc[omics_data.index.str.contains("LL")]
            template_model.exchanges.EX_C00205__dra.bounds = (-150*2.99, -150*2.99)
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


def get_best_thresholds(df: pd.DataFrame):
    dist = rank_models_difference(df)
    sorted_thresholds = sorted(dist, key=lambda k: dist[k], reverse=True)
    return sorted_thresholds[0]


def reconstruct_models(config_path: str):
    params = json.load(open(rf"{config_path}", "r"))
    media_path = params.get("MEDIA_PATH")
    media_sheet = params.get("MEDIA_SHEET")
    model_path = params.get("MODEL_PATH")
    consistent_model_path = params.get("CONSISTENT_MODEL_PATH")
    media = pd.read_excel(media_path, index_col=0, sheet_name=None, engine='openpyxl')[media_sheet].to_dict(orient='index')
    medium_conditions = {'f2 medium': {key: (value['LB'], value['UB']) for key, value in media.items()}}
    params = json.load(open(rf"{config_path}", "r"))
    template_model = load_model(model_path=model_path, consistent_model_path=consistent_model_path, medium_conditions=medium_conditions)
    # reaction_id  = "PRISM_solar_litho__extr"
    reaction_id = 'PRISM_white_LED__extr'
    for reaction in template_model.reactions:
        if reaction.id.startswith("PRISM") and reaction.id != reaction_id:
            reaction.bounds = (0, 0)
    template_model.reactions.get_by_id(reaction_id).bounds = (0, 10000)
    for algorithm in params.get("ALGORITHMS"):
        file_path = f'{algorithm}_{params["THRESHOLDING_STRATEGY"].replace(" ", "_")}_' \
                    f'{params["GLOBAL_THRESHOLD_UPPER"]}_{params["GLOBAL_THRESHOLD_LOWER"]}_{params["LOCAL_THRESHOLD"]}.csv'
        df = pd.read_csv(os.path.join(params.get("TROPPO_RESULTS_PATH"), 'integration',file_path), index_col=0)
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


if __name__ == "__main__":
    # config = "nextflow/PRJNA589063.json"
    config = "nextflow/PRJNA495151.json"
    integrate(config)
    reconstruct_models(config)