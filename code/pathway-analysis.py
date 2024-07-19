import json

import pandas as pd

from gsmmutils import MyModel, DATA_PATH
from gsmmutils.api.kegg import search_pathway_map_id
import webbrowser


def load_data(filenames):
    """
    Load data from files
    :param filenames:
    :return:
    """
    res = {}
    for filename in filenames:
        res[filename.split("\\")[-1].split(".csv")[0]] = [pd.read_csv(filename, index_col=0), pd.read_csv(filename.replace("_results", ""))]
    return res


def paint_KEGG_pathway(reaction_ids, pathway_name, model):
    pathway_id = search_pathway_map_id(pathway_name)
    # url = f"https://www.kegg.jp/pathway/{pathway_id}"
    if pathway_id:
        url = f"https://www.kegg.jp/kegg-bin/show_pathway?map={pathway_id}&multi_query="
        for reaction, info in reaction_ids.iterrows():
            reaction_id = None
            if not reaction.startswith("R"):
                if "ec-code" in model.reactions.get_by_id(reaction).annotation:
                    reaction_id = model.reactions.get_by_id(reaction).annotation["ec-code"]
                    if type(reaction_id) == list:
                        if info['FC'] > 0:
                            reaction_id = color_from_fc(info['FC'], '%20green%0A'.join(reaction_id))
                        else:
                            reaction_id = color_from_fc(info['FC'], '%20red%0A"'.join(reaction_id))
                    else:
                        reaction_id = color_from_fc(info['FC'], reaction_id)
                else:
                    print(f"Reaction {reaction} doesn't have ec-code annotation")
            else:
                reaction_id = color_from_fc(info['FC'], reaction.split("__")[0])
            if reaction_id:
                url += reaction_id
        print(url)
        webbrowser.open(url)

def color_from_fc(fc, reaction_id):
    if fc > 0:
        return reaction_id + "%20green%0A"
    else:
        return reaction_id + "%20red%0A"


def get_reactions_pathway_map(model_path: str):
    model = MyModel(model_path, "e_Biomass__cytop")
    json.dump(model.reactions_pathways_map, open(rf"{DATA_PATH}/omics\reactions_pathways_map.json", "w"))
    json.dump(model.pathway_reactions_map, open(rf"{DATA_PATH}/omics\pathways_reactions_map.json", "w"))
    reactions_pathways_map = json.load(open(rf"{DATA_PATH}/omics\reactions_pathways_map.json", "r"))
    pathways_reactions_map = json.load(open(rf"{DATA_PATH}/omics\pathways_reactions_map.json", "r"))
    return reactions_pathways_map, pathways_reactions_map, model


def get_der_by_pathway(pathway_name, pathway_map, data, model):
    reactions = data.loc[data.index.isin(pathway_map[pathway_name])]
    paint_KEGG_pathway(reactions, pathway_name, model)


if __name__ == '__main__':
    data = load_data([
        r"C:\Users\Bisbii\PythonProjects\gsmmutils\data\omics\nacl_fastcore_results.csv",
        r"C:\Users\Bisbii\PythonProjects\gsmmutils\data\omics\nacl_gimme_results.csv",
        #r"C:\Users\Bisbii\PythonProjects\gsmmutils\data\omics\h2o2_fastcore_results.csv",
        #r"C:\Users\Bisbii\PythonProjects\gsmmutils\data\omics\h2o2_gimme_results.csv",
        # r"C:\Users\Bisbii\PythonProjects\gsmmutils\data\omics\ml_ll_results.csv",
        # r"C:\Users\Bisbii\PythonProjects\gsmmutils\data\omics\hl_ml_results.csv",
        # r"C:\Users\Bisbii\PythonProjects\gsmmutils\data\omics\hl_ll_results.csv"
    ])
    reactions_pathways_map, pathways_reactions_map, model = get_reactions_pathway_map("../data/ngaditana/models/model_ng.xml")
    for key in data.keys():
        for pathway in data[key][1].Pathways.tolist():
            if pathway != "Fatty acid metabolism":
                get_der_by_pathway(pathway, pathways_reactions_map, data[key][0], model)
