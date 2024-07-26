import json
import os.path
import webbrowser
from os.path import join

import gseapy as gp
import pandas as pd
from gsmmutils import MyModel, DATA_PATH
from gsmmutils.api.kegg import search_pathway_map_id


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
                            reaction_id = color_from_fc(info['FC'], '%20red%0A'.join(reaction_id))
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
    else:
        print(f"Pathway not found: {pathway_name}")

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


def reaction_analysis():
    data = load_data([
        join(RESULTS_PATH, "dfa/C5_fastcore_P5_fastcore_results.csv")
    ])
    reactions_pathways_map, pathways_reactions_map, model = get_reactions_pathway_map(join(DATA_PATH, "models/model_ng.xml"))
    for key in data.keys():
        for pathway in data[key][1].Pathways.tolist():
            if pathway != "Fatty acid metabolism":
                get_der_by_pathway(pathway, pathways_reactions_map, data[key][0], model)


def pathway_enrichment(gene_symbols, custom_mapping):
    with open(join(RESULTS_PATH, 'deg/custom_pathway.gmt'), 'w') as gmt_file:
        for pathway, group in custom_mapping.groupby('PathwayID'):
            genes = '\t'.join(group['GeneID'])
            gmt_file.write(f"{pathway}\tCustom Pathway\t{genes}\n")

    enrichr_results = gp.enrichr(
        gene_list=gene_symbols,
        gene_sets=join(RESULTS_PATH, 'deg/custom_pathway.gmt'),
        organism='hs',
        outdir=join(RESULTS_PATH, 'deg/enrichr'),
        cutoff=0.05
    )
    enrichr_results.results.to_csv('custom_pathway_enrichment.csv', index=False)

def get_gene_pathway_map(model):
    gene_pathway_ids = {}
    pathway_ids = {}
    pathway_id_c = 0
    if not os.path.exists(join(RESULTS_PATH, "pathway_ids.json")):
        for pathway in model.groups:
            kegg_id = search_pathway_map_id(pathway.name)
            if kegg_id:
                pathway_ids[pathway.name] = kegg_id
            else:
                pathway_ids[pathway.name] = f"ng{pathway_id_c}"
                pathway_id_c += 1
        with open(join(RESULTS_PATH, "pathway_ids.json"), "w") as f:
            json.dump(pathway_ids, f, indent=4)
    else:
        pathway_ids = json.load(open(join(RESULTS_PATH, "pathway_ids.json"), "r"))
    for gene in model.genes:
        gene_pathway_ids[gene.id] = []
        for pathway in model.genes_pathways_map[gene.id]:
            gene_pathway_ids[gene.id].append(pathway_ids[pathway])

    gene_pathway_ids = pd.DataFrame.from_dict(gene_pathway_ids, orient='index')
    gene_pathway_ids["GeneID"] = gene_pathway_ids.index
    gene_pathway_ids.reset_index(drop=True, inplace=True)
    gene_pathway_ids = gene_pathway_ids.melt(id_vars="GeneID", value_name = "PathwayID")
    gene_pathway_ids = gene_pathway_ids.dropna(subset=['PathwayID'])
    gene_pathway_ids.to_csv(join(RESULTS_PATH, "deg/gene_pathway.csv"), sep="\t")
    return gene_pathway_ids


def gene_analysis():
    data =  pd.read_csv(join(RESULTS_PATH, "deg/degs.tsv"), sep="\t", index_col=0)
    data.index = [e.replace("-", "_") for e in data.index]
    _, _, model = get_reactions_pathway_map(join(DATA_PATH, "models/model_ng.xml"))
    model.get_genes_pathways_map()
    if not os.path.exists(join(RESULTS_PATH, "deg/gene_pathway.csv")):
        gene_pathway_ids = get_gene_pathway_map(data, model)
    else:
        gene_pathway_ids = pd.read_csv(join(RESULTS_PATH, "deg/gene_pathway.csv"), sep="\t", index_col=0)
    pathway_enrichment(data.index.to_list(), gene_pathway_ids)




if __name__ == '__main__':
    DATA_PATH = "/home/ecunha/omics-integration/data/ngaditana"
    RESULTS_PATH = "/home/ecunha/omics-integration/results/ngaditana/PRJNA589063"
    reaction_analysis()
    # gene_analysis()
