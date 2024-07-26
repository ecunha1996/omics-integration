import os

import pandas as pd

from gsmmutils.graphics.plot import clustermap
from gsmmutils.io import read_csv
from gsmmutils.model.COBRAmodel import MyModel
import pickle

pd.set_option('display.max_rows', 500)
pd.set_option('display.max_columns', 500)


def get_omics_object():
    omics_light = pickle.load(open("light/omics_light.pkl", "rb"))
    omics_nacl_h2o2_sorb = pickle.load(open("nacl_h2o2_sorb/nacl_h2o2_sorb.pkl", "rb"))
    return omics_light, omics_nacl_h2o2_sorb


def gimme(results_file, omics, model, method="GIMME"):
    omics.get_flux_change(method_1="GIMME", condition_1="control", condition_2="nacl", method_2="GIMME", threshold=0.67)
    omics.get_flux_change(method_1="GIMME", condition_1="control", condition_2="sorb", method_2="GIMME", threshold=0.67)
    omics.get_flux_change(method_1="GIMME", condition_1="control", condition_2="h2o2", method_2="GIMME", threshold=0.67)
    df = pd.concat([value.dataframe.rename(columns={'Flux rate': key}) for key, value in omics.specific_models['GIMME'].items()], axis=1)
    to_keep = list(set(list(omics.flux_change['GIMME_control_GIMME_nacl'].keys()) + list(omics.flux_change['GIMME_control_GIMME_sorb'].keys()) +
                       list(omics.flux_change['GIMME_control_GIMME_h2o2'].keys())))
    df = df.loc[to_keep]
    clustermap(df, title="GIMME")
    get_reaction_flux_expression_data(model, omics, method)


def get_reaction_flux_expression_data(model, omics, method):
    for key, value in omics.flux_change.items():
        res = {}
        for reaction, data in value.items():
            if 'ec-code' in model.reactions.get_by_id(reaction).annotation.keys():
                if type(model.reactions.get_by_id(reaction).annotation['ec-code']) is list:
                    ecnumber = '; '.join(model.reactions.get_by_id(reaction).annotation['ec-code'])
                else:
                    ecnumber = model.reactions.get_by_id(reaction).annotation['ec-code']
            else:
                ecnumber = ''
            pathway = '; '.join(model.reactions_pathways_map[reaction])
            res[reaction] = [data, ecnumber, pathway]
        as_df = pd.DataFrame.from_dict(res, orient='index', columns=['Flux rate', 'EC number', 'Pathway'])
        as_df.to_csv(f"nacl_h2o2_sorb/{key}.csv")
    pathway_counts_nacl_gimme = count_reaction_by_pathway(omics.flux_change[f'{method}_control_{method}_nacl'], model)
    print(dict(sorted(pathway_counts_nacl_gimme.items(), key=lambda item: item[1], reverse=True)))


def count_reaction_by_pathway(flux_changes: dict, model: MyModel):
    pathway_counts = {}
    for reaction in flux_changes.keys():
        for pathway in model.reactions_pathways_map[reaction]:
            if pathway in pathway_counts.keys():
                pathway_counts[pathway] += 1
            else:
                pathway_counts[pathway] = 1
    return pathway_counts


def imat():
    pass


def fastcore(model):
    data = read_csv(r"nacl_results.csv")
    print(data.head())


def tinit():
    print(os.getcwd())
    model = MyModel("nacl_h2o2_sorb\sorb\Dsalina_sorb_Local2_2_4_4_tinit_t0.xml", "e_Biomass__cytop")


def sampling_analysis(model, results_file, pathway=None):
    data = read_csv(results_file)
    if pathway:
        model.get_pathway_reactions_map()
        reactions_in_pathway = model.pathway_reactions_map[pathway]
        data_pathway = data.loc[data['Reaction'].isin(reactions_in_pathway)]
        for reaction in data_pathway['Reaction']:
            if 'ec-code' in model.reactions.get_by_id(reaction).annotation.keys():
                if type(model.reactions.get_by_id(reaction).annotation['ec-code']) is list:
                    ecnumber = '; '.join(model.reactions.get_by_id(reaction).annotation['ec-code'])
                else:
                    ecnumber = model.reactions.get_by_id(reaction).annotation['ec-code']
            else:
                ecnumber = ''
            pathway = '; '.join(model.reactions_pathways_map[reaction])
            data_pathway.loc[data_pathway['Reaction'] == reaction, 'EC number'] = ecnumber
            data_pathway.loc[data_pathway['Reaction'] == reaction, 'Pathway'] = pathway
        under_expressed = data_pathway.loc[data_pathway['FC'] < 0]
        over_expressed = data_pathway.loc[data_pathway['FC'] > 0]
        print(under_expressed[['Reaction', 'EC number']].sort_values(by=['EC number'], ascending=True))
        print("#" * 100)
        print(over_expressed[['Reaction', 'EC number']].sort_values(by=['EC number'], ascending=True))


if __name__ == '__main__':
    os.chdir("../../data/omics")
    # omics_light, omics_nacl_h2o2_sorb = get_omics_object()
    model = MyModel("../models/model_with_trials.xml", "e_Biomass__cytop")
    # gimme("context_specific_models.xlsx", omics_nacl_h2o2_sorb, model)
    # fastcore(model)
    # tinit()
    sampling_analysis(model, 'nacl_results.csv', 'Pyrimidine metabolism')
