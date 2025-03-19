

import json
import os.path
import webbrowser
from collections import OrderedDict
from copy import deepcopy
from os.path import join

import pandas as pd
from gsmmutils.model import MyModel
from gsmmutils.api.kegg import search_pathway_map_id
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import sys
sys.path.append("../../")
from utils import load_samples
plt.rcParams["font.family"] = "Arial"


def load_data(filenames):
    """
    Load data from files
    :param filenames:
    :return:
    """
    res = {}
    for key, filename in filenames.items():
        tmp = pd.read_csv(filename, index_col=0, sep="\t")
        as_df = pd.read_csv(filename.replace("_all_results.tsv", ".csv")).sort_values(by="P-value_adj", ascending=True) #.iloc[:20, ]
        res[key] = [(tmp.loc[(abs(tmp['FC']) >= 0.0) & (tmp['Padj'] <= 0.05), :]), as_df]
    # rename "Ubiquinone and other terpenoid-quinone biosynthesis" to "Quinone biosynthesis"
    for key, df in res.items():
        df[1].Pathways = df[1].Pathways.replace("Ubiquinone and other terpenoid-quinone biosynthesis", "Quinone biosynthesis")
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
        print(f"Pathway found: {pathway_name}")
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
    json.dump(model.reactions_pathways_map, open(rf"{DATA_PATH}/reactions_pathways_map.json", "w"))
    json.dump(model.pathway_reactions_map, open(rf"{DATA_PATH}/pathways_reactions_map.json", "w"))
    reactions_pathways_map = json.load(open(rf"{DATA_PATH}/reactions_pathways_map.json", "r"))
    pathways_reactions_map = json.load(open(rf"{DATA_PATH}/pathways_reactions_map.json", "r"))
    # rename "Ubiquinone and other terpenoid-quinone biosynthesis" to "Quinone biosynthesis"
    pathways_reactions_map["Quinone biosynthesis"] = pathways_reactions_map.pop("Ubiquinone and other terpenoid-quinone biosynthesis")
    for reaction, pathways in reactions_pathways_map.items():
        if "Ubiquinone and other terpenoid-quinone biosynthesis" in pathways:
            pathways.remove("Ubiquinone and other terpenoid-quinone biosynthesis")
            pathways.append("Quinone biosynthesis")
    return reactions_pathways_map, pathways_reactions_map, model


def get_der_by_pathway(pathway_name, pathway_map, data, model):
    reactions = data.loc[data.index.isin(pathway_map[pathway_name])]
    paint_KEGG_pathway(reactions, pathway_name, model)


def reaction_analysis():
    data = load_data({
        # join(RESULTS_PATH, "dfa/C5_fastcore_0.4_loopless_ACHR_P5_fastcore_0.4_loopless_ACHR_results.csv")
        # join(RESULTS_PATH, "dfa/control_fastcore_0.6_loopless_ACHR_nacl_fastcore_0.6_loopless_ACHR_all_results.tsv")
        # join(RESULTS_PATH, "dfa/control_gimme_0.4_loopless_ACHR_nacl_gimme_0.4_loopless_ACHR_all_results.tsv")
        # "Control vs H2O2 GIMME": join(RESULTS_PATH, "dfa/control_gimme_0.4_ACHR_h2o2_gimme_0.4_ACHR_all_results.tsv"),
        # "Control vs NaCl GIMME": join(RESULTS_PATH, "dfa/control_gimme_0.4_ACHR_nacl_gimme_0.4_ACHR_all_results.tsv"),
        # "Control vs Sorb GIMME": join(RESULTS_PATH, "dfa/control_gimme_0.4_ACHR_sorb_gimme_0.4_ACHR_all_results.tsv"),

        # "Control vs H2O2 FASTCORE": join(RESULTS_PATH, "dfa/control_fastcore_0.6_ACHR_h2o2_fastcore_0.6_ACHR_all_results.tsv"),
        # "Control vs NaCl FASTCORE": join(RESULTS_PATH, "dfa/control_fastcore_0.6_ACHR_nacl_fastcore_0.6_ACHR_all_results.tsv"),
        # #
        # "LL vs ML FASTCORE": join(RESULTS_PATH, "dfa/LL_fastcore_0.65_ACHR_ML_fastcore_0.65_ACHR_all_results.tsv"),
        "HL vs LL FASTCORE": join(RESULTS_PATH, "dfa/LL_fastcore_0.65_ACHR_HL_fastcore_0.65_ACHR_all_results.tsv"),

        # "LL vs ML FASTCORE": join(RESULTS_PATH, "dfa/LL_fastcore_0.65_ACHR_ML_fastcore_0.65_ACHR_all_results.tsv"),
        # "HL vs LL FASTCORE": join(RESULTS_PATH, "dfa/LL_fastcore_0.65_ACHR_samples_HL_fastcore_0.65_ACHR_samples_all_results.tsv"),

        # "Control vs H2O2 FASTCORE": join(RESULTS_PATH, "../PRJNA437866/dfa/control_fastcore_0.6_ACHR_h2o2_fastcore_0.6_ACHR_all_results.tsv"),
        # "Control vs NaCl FASTCORE": join(RESULTS_PATH, "../PRJNA437866/dfa/control_fastcore_0.6_ACHR_nacl_fastcore_0.6_ACHR_all_results.tsv"),
        # "Control vs SORB FASTCORE": join(RESULTS_PATH, "../PRJNA437866/dfa/control_fastcore_0.6_ACHR_sorb_fastcore_0.6_ACHR_all_results.tsv"),


        # "Control vs H2O2 FASTCORE": join(RESULTS_PATH, "../PRJNA437866/dfa/control_fastcore_0.8_ACHR_h2o2_fastcore_0.8_ACHR_all_results.tsv"),

        # "LL vs ML FASTCORE": join(RESULTS_PATH, "dfa/LL_fastcore_0.6_ACHR_ML_fastcore_0.6_ACHR_all_results.tsv"),
        # "HL vs LL FASTCORE": join(RESULTS_PATH, "dfa/LL_fastcore_0.6_ACHR_HL_fastcore_0.6_ACHR_all_results.tsv"),

        # "LL vs ML GIMME": join(RESULTS_PATH, "dfa/LL_gimme_0.4_ACHR_ML_gimme_0.4_ACHR_all_results.tsv"),
         "LL vs HL GIMME": join(RESULTS_PATH, "dfa/LL_gimme_0.4_ACHR_HL_gimme_0.4_ACHR_all_results.tsv"),
    })
    reactions_pathways_map, pathways_reactions_map, model = get_reactions_pathway_map(join(DATA_PATH, "models/model_ds.xml"))

    common_pathways = set() #set(data[list(data.keys())[0]][1].Pathways.tolist()).intersection(set(data[list(data.keys())[1]][1].Pathways.tolist())) #
    fig, ax = plt.subplots(1, 2, figsize=(12, 8))
    for i, (condition, df) in enumerate(data.items()):
        dfa_paths = data[list(data.keys())[i]][1].Pathways.tolist()
        ordered_pathway_map = OrderedDict()
        for pathway in dfa_paths:
            if pathway in pathways_reactions_map:
                ordered_pathway_map[pathway] = pathways_reactions_map[pathway]
        ordered_pathway_map = OrderedDict(reversed(list(ordered_pathway_map.items())))
        fc_plot(ax[i], data[list(data.keys())[i]][0], ordered_pathway_map,
                condition, common_pathways, i)

    ax[0].text(-0.4, 1.05, 'A', fontsize=10, transform=ax[0].transAxes, fontweight='bold', va='top', ha='right')
    ax[1].text(-0.4, 1.05, 'B', fontsize=10, transform=ax[1].transAxes, fontweight='bold', va='top', ha='right')

    plt.tight_layout()
    plt.savefig(f"{RESULTS_PATH}/dsalina_light.pdf", bbox_inches='tight', format='pdf', dpi=1200)
    plt.show()
    plt.close()

    #
    # get_der_by_pathway("Carbon fixation in photosynthetic organisms", pathways_reactions_map, data[list(data.keys())[0]][0], model)
    # get_der_by_pathway("Carbon fixation in photosynthetic organisms", pathways_reactions_map, data[list(data.keys())[0]][1], model)
    # get_der_by_pathway("Glycine, serine and threonine metabolism", pathways_reactions_map, data[list(data.keys())[0]][0], model)
    # get_der_by_pathway("Pentose phosphate pathway", pathways_reactions_map, data[list(data.keys())[0]][0], model)
    get_der_by_pathway("Citrate cycle (TCA cycle)", pathways_reactions_map, data[list(data.keys())[0]][0], model)
    #get_der_by_pathway("Cysteine and methionine metabolism", pathways_reactions_map, data[list(data.keys())[0]][0], model)
    # get_der_by_pathway("Pyruvate metabolism", pathways_reactions_map, data[list(data.keys())[0]][0], model)
    # get_der_by_pathway("Glycolysis / Gluconeogenesis", pathways_reactions_map, data[list(data.keys())[0]][0], model)
    # get_der_by_pathway("Carotenoid biosynthesis", pathways_reactions_map, data[list(data.keys())[0]][0], model)
    # get_der_by_pathway("Carotenoid biosynthesis", pathways_reactions_map, data[list(data.keys())[1]][0], model)
    # get_der_by_pathway("Starch and sucrose metabolism", pathways_reactions_map, data[list(data.keys())[0]][0], model)
    # get_der_by_pathway("Fatty acid degradation", pathways_reactions_map, data[list(data.keys())[0]][0], model)
    get_der_by_pathway("Glyoxylate and dicarboxylate metabolism", pathways_reactions_map, data[list(data.keys())[0]][0], model)
    # get_der_by_pathway("Carbon fixation pathways in prokaryotes", pathways_reactions_map, data[list(data.keys())[0]][0], model)
    # get_der_by_pathway("Valine, leucine and isoleucine biosynthesis", pathways_reactions_map, data[list(data.keys())[0]][0], model)
    # for key in data.keys():
    #     for pathway in data[key][1].Pathways.tolist():
    #         if pathway != "Fatty acid metabolism":
    #             get_der_by_pathway(pathway, pathways_reactions_map, data[key][0], model)


def reaction_analysis_all():
    data = load_data({

        "LL vs ML FASTCORE": join(RESULTS_PATH, "dfa/LL_fastcore_0.65_ACHR_ML_fastcore_0.65_ACHR_all_results.tsv"),
        "LL vs ML GIMME": join(RESULTS_PATH, "dfa/LL_gimme_0.4_ACHR_ML_gimme_0.4_ACHR_all_results.tsv"),
        "HL vs LL FASTCORE": join(RESULTS_PATH, "dfa/LL_fastcore_0.65_ACHR_HL_fastcore_0.65_ACHR_all_results.tsv"),
        "HL vs LL GIMME": join(RESULTS_PATH, "dfa/LL_gimme_0.4_ACHR_HL_gimme_0.4_ACHR_all_results.tsv"),

        # "Control vs H2O2 FASTCORE": join(RESULTS_PATH, "../PRJNA437866/dfa/control_fastcore_0.6_ACHR_h2o2_fastcore_0.6_ACHR_all_results.tsv"),
        # "Control vs NaCl FASTCORE": join(RESULTS_PATH, "../PRJNA437866/dfa/control_fastcore_0.6_ACHR_nacl_fastcore_0.6_ACHR_all_results.tsv"),

    })
    reactions_pathways_map, pathways_reactions_map, model = get_reactions_pathway_map(join(DATA_PATH, "models/model_ds.xml"))

    common_pathways = set() #set(data[list(data.keys())[0]][1].Pathways.tolist()).intersection(set(data[list(data.keys())[1]][1].Pathways.tolist())) #
    fig, axs = plt.subplots(2, 2, figsize=(7.08, 9))
    for i, (condition, df) in enumerate(data.items()):
        if i % 2 == 0:
            common_pathways = set(data[list(data.keys())[i]][1].Pathways.tolist()).intersection(set(data[list(data.keys())[i+1]][1].Pathways.tolist()))
        dfa_paths = data[list(data.keys())[i]][1].Pathways.tolist()
        ordered_pathway_map = OrderedDict()
        for pathway in dfa_paths:
            if pathway in pathways_reactions_map:
                ordered_pathway_map[pathway] = pathways_reactions_map[pathway]
        ordered_pathway_map = OrderedDict(reversed(list(ordered_pathway_map.items())))
        fc_plot(axs[i//2, i%2], data[list(data.keys())[i]][0], ordered_pathway_map,
                condition, common_pathways, i, remove_first=False)

    axs[0][0].text(-0.4, 1.05, 'A', fontsize=10, transform=axs[0][0].transAxes, fontweight='bold', va='top', ha='right')
    axs[0][1].text(-0.4, 1.05, 'B', fontsize=10, transform=axs[0][1].transAxes, fontweight='bold', va='top', ha='right')
    axs[1][0].text(-0.4, 1.05, 'C', fontsize=10, transform=axs[1][0].transAxes, fontweight='bold', va='top', ha='right')
    axs[1][1].text(-0.4, 1.05, 'D', fontsize=10, transform=axs[1][1].transAxes, fontweight='bold', va='top', ha='right')
    plt.tight_layout()
    fig.subplots_adjust(hspace=0.12)
    plt.savefig(f"{RESULTS_PATH}/dsalina_light_all.pdf", bbox_inches='tight', format='pdf', dpi=1200)
    plt.show()
    plt.close()





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
    # if not os.path.exists(join(RESULTS_PATH, "pathway_ids.json")):
    for pathway in model.groups:
        kegg_id =  pathway.name.replace(" ", "_") #search_pathway_map_id(pathway.name)
        if kegg_id:
            pathway_ids[pathway.name] = kegg_id
        else:
            pathway_ids[pathway.name] = f"ng{pathway_id_c}"
            pathway_id_c += 1
    with open(join(RESULTS_PATH, "pathway_ids.json"), "w") as f:
        json.dump(pathway_ids, f, indent=4)
    # else:
    #     pathway_ids = json.load(open(join(RESULTS_PATH, "pathway_ids.json"), "r"))
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
    _, _, model = get_reactions_pathway_map(join(DATA_PATH, "models/model_ds.xml"))
    model.get_genes_pathways_map()
    use_cache = False
    if not os.path.exists(join(RESULTS_PATH, "deg/gene_pathway.csv")) or not use_cache:
        gene_pathway_ids = get_gene_pathway_map(model)
    else:
        gene_pathway_ids = pd.read_csv(join(RESULTS_PATH, "deg/gene_pathway.csv"), sep="\t", index_col=0)
    pathway_enrichment(data.index.to_list(), gene_pathway_ids)

def fc_plot(ax, data, pathways_map, condition_name, common_pathways, plot_number, remove_first=True):

    np.random.seed(42)
    values = {}

    for pathway, reactions in pathways_map.items():
        values[pathway] = data.loc[data.index.isin(reactions)].FC.tolist()
    bar_width = 0.4

    for i, (pathway, value) in enumerate(values.items()):
        positive_vals = len([v for v in value if v > 0])
        negative_vals = len([v for v in value if v < 0])
        ax.barh(y=pathway, width=negative_vals, height=bar_width, color='r', label='Negative FC' if i == 0 else "")
        ax.barh(y=pathway, width=positive_vals, height=bar_width, color='g', label='Positive FC' if i == 0 else "", left = negative_vals)

    ax.set_yticklabels([split_label(label, 50) for label in values.keys()])
    for lab in ax.get_yticklabels():
        if lab.get_text().replace("\n", "") in common_pathways or lab.get_text().replace("\n", " ") in common_pathways:
            lab.set_fontweight('bold')
        lab.set_fontsize(7)

    if remove_first and plot_number == 0:
        ax.set_xlabel('')
        ax.set_xticks([])
        ax.set_xticklabels([])
        ax.legend(fontsize=7)
    else:
        ax.set_xlabel('Number of Reactions', fontsize=7)
        ax.tick_params(axis='x', labelsize=7)
        ax.legend().remove()

# def fc_plot(fig, ax, data, pathways_map, condition_name, common_pathways, plot_number):
#
#     np.random.seed(42)
#     values = {}
#
#     for pathway, reactions in pathways_map.items():
#         values[pathway] = data.loc[data.index.isin(reactions)].FC.tolist()
#     bar_width = 0.6
#
#     # plt.figure(figsize=(6.4, 6.5))
#     # for i, (pathway, value) in enumerate(values.items()):
#     #     plt.scatter(value, [pathway] * len(value), label=pathway, alpha=0.6)
#     #
#     # plt.xlabel('FC')
#     # plt.title(condition_name)
#     # plt.grid(True)
#     # plt.yticks(ha='right')
#     # plt.tight_layout()
#     # plt.savefig(f"/home/ecunha/omics-integration/results/ngaditana/PRJNA589063/{condition_name}_scatter.png", bbox_inches='tight', dpi=800)
#     # plt.show()
#     # plt.close()
#
#     # plt.figure(figsize=(6.4, 4.8))
#     # for i, (pathway, value) in enumerate(values.items()):
#     #     positive_vals = len([v for v in value if v > 0])
#     #     negative_vals = len([v for v in value if v < 0])
#     #     plt.barh(y=i - bar_width, width=negative_vals, height=bar_width, color='r', label='Negative')
#     #     plt.barh(y=pathway, width=positive_vals, height=bar_width, color='g', label='Positive')
#     #
#     # plt.legend()
#     # plt.xlabel('Number of Reactions')
#     # plt.title(condition_name)
#     # plt.grid(True)
#     # plt.yticks(ha='right')
#     # plt.tight_layout()
#     # plt.savefig(f"/home/ecunha/omics-integration/results/ngaditana/PRJNA589063/{condition_name}_bar.png", bbox_inches='tight', dpi=600)
#     # plt.show()
#     # plt.close()
#
#     for i, (pathway, value) in enumerate(values.items()):
#         positive_vals = len([v for v in value if v > 0])
#         negative_vals = len([v for v in value if v < 0])
#         ax.barh(y=pathway, width=negative_vals, height=bar_width, color='r', label='Negative FC' if i == 0 else "")
#         ax.barh(y=pathway, width=positive_vals, height=bar_width, color='g', label='Positive FC' if i == 0 else "", left = negative_vals)
#
#     ax.set_yticklabels([split_label(label, 32) for label in values.keys()], linespacing=0.85)
#
#     for lab in ax.get_yticklabels():
#         if lab.get_text().replace("\n", "") in common_pathways or lab.get_text().replace("\n", " ") in common_pathways:
#             lab.set_fontweight('bold')
#         lab.set_fontsize(7)
#
#     ax.legend(fontsize=7)
#     ax.set_xlabel('Number of Reactions', fontsize=7)
#     # ax.set_title(' '.join(condition_name.split()[0:-1]))
#     # add a letter on the topleft of the plot
#     # ax.grid(True)
#     ax.tick_params(axis='x', labelsize=7)
#     number_letter_map = {0: 'A', 1: 'B'}
#     ax.text(-0.5, 1.04, number_letter_map[plot_number], fontsize=10, transform=ax.transAxes, fontweight='bold', va='top', ha='right')
#     # ax.set_yticks(ha='right')

def split_labels(labels, max_line_length=45):
    split_labels = []
    for label in labels:
        words = label.split()
        new_label = ""
        current_line = ""
        for word in words:
            if len(current_line) + len(word) + 1 > max_line_length:
                new_label += current_line + "\n"
                current_line = word
            else:
                if current_line:
                    current_line += " " + word
                else:
                    current_line = word
        new_label += current_line
        split_labels.append(new_label)
    return split_labels

def split_label(label, max_line_length=25):
    words = label.split()
    new_label = ""
    current_line = ""
    for word in words:
        if len(current_line) + len(word) + 1 > max_line_length:
            new_label += current_line + "\n"
            current_line = word
        else:
            if current_line:
                current_line += " " + word
            else:
                current_line = word
    new_label += current_line
    return new_label

def get_enzymes_map(model, bmgr):
    labels = []
    ec_number_enzyme_map = {"3.1.1.3": "LIP", "2.3.1.20": "DGAT", "2.3.1.158": "PDAT", "1.14.19.25": "FAD3",
                            "3.1.1.23": "MLIP", "3.1.3.36": "INPP5K", "2.7.1.68": "PIP5K", "2.7.8.11": "PIS", "3.1.3.4": "PAH1",
                            "1.14.19.30": "D5Des", "2.3.1.22": "LPGAT1", "1.14.19.47": "DES6", "1.14.19.22": "FAD2",
                            "1.14.19.25,1.14.19.35": "FAD3", "2.7.1.107": "DGKG", '1.14.19.25,1.14.19.36': "FAD3",
                            "1.14.19.45": "desA", "2.3.1.15": "GPAT", "2.3.1.51": "LPAT"}
    for reaction in bmgr.index.tolist():
        ec_number = model.reactions.get_by_id(reaction).annotation.get("ec-code", "")
        if isinstance(ec_number, list):
            ec_number = ','.join(ec_number)
        labels.append(ec_number_enzyme_map.get(ec_number, ""))
    return labels

def plot_heatmap():
    data = load_samples(RESULTS_PATH, samples=["C3", "C5", "N3", "N5", "P3", "P5"])
    g = sns.clustermap(data['all'], cmap="viridis", yticklabels=False)
    plt.setp(g.ax_heatmap.get_xticklabels(), rotation=45, fontsize=14)
    plt.savefig(join(RESULTS_PATH, "ngaditana_heatmap_complete.png"), dpi=600)
    gimme = data['gimme']
    model = MyModel(f"{DATA_PATH}/models/model_ng.xml", "e_Biomass__cytop")
    bmgr = gimme.loc[(gimme.index.str.startswith("BMGR")) & (
            (gimme.index.isin(model.pathway_reactions_map["BOIMMG (TAG)"])) |
            ((gimme.index.isin(model.pathway_reactions_map["BOIMMG (DAG)"]))) |
            ((gimme.index.isin(model.pathway_reactions_map["BOIMMG (PA)"])))
    )]
    bmgr = bmgr.filter(regex="N3|C3|P3|C5|P5|N5_gimme")
    bmgr = bmgr.loc[~bmgr.index.str.contains("mem")]
    row_std = bmgr.std(axis=1)
    bmgr = bmgr[row_std != 0]
    labels = get_enzymes_map(model, bmgr)
    tmp = deepcopy(bmgr)
    tmp.index = labels
    grouped_df = tmp.groupby(tmp.index).mean()
    g = sns.clustermap(grouped_df, cmap="viridis", z_score=0)
    plt.savefig(join(RESULTS_PATH, "dsalina_heatmap_bmgr.png"), dpi=600)


def reaction_analysis_stress():
    data = load_data({
        "Control vs NaCl FASTCORE": join(RESULTS_PATH, "../PRJNA437866/dfa/control_fastcore_0.6_ACHR_nacl_fastcore_0.6_ACHR_all_results.tsv"),
        "Control vs NaCl GIMME": join(RESULTS_PATH, "dfa/control_gimme_0.4_ACHR_nacl_gimme_0.4_ACHR_all_results.tsv"),
        "Control vs Sorb FASTCORE": join(RESULTS_PATH, "../PRJNA437866/dfa/control_fastcore_0.6_ACHR_sorb_fastcore_0.6_ACHR_all_results.tsv"),
        "Control vs Sorb GIMME": join(RESULTS_PATH, "dfa/control_gimme_0.4_ACHR_sorb_gimme_0.4_ACHR_all_results.tsv"),
        "Control vs H2O2 FASTCORE": join(RESULTS_PATH, "../PRJNA437866/dfa/control_fastcore_0.6_ACHR_h2o2_fastcore_0.6_ACHR_all_results.tsv"),
        "Control vs H2O2 GIMME": join(RESULTS_PATH, "dfa/control_gimme_0.4_ACHR_h2o2_gimme_0.4_ACHR_all_results.tsv"),

    })
    reactions_pathways_map, pathways_reactions_map, model = get_reactions_pathway_map(join(DATA_PATH, "models/model_ds.xml"))
    common_pathways = set() #set(data[list(data.keys())[0]][1].Pathways.tolist()).intersection(set(data[list(data.keys())[1]][1].Pathways.tolist())) #
    fig, ax = plt.subplots(3, 2, figsize=(7.08, 8))
    j = 0
    for i, (condition, df) in enumerate(data.items()):
        if i % 2 == 0:
            common_pathways = set(data[list(data.keys())[i]][1].Pathways.tolist()).intersection(set(data[list(data.keys())[i+1]][1].Pathways.tolist()))
        if j == 2:
            j = 0
        dfa_paths = data[list(data.keys())[i]][1].Pathways.tolist()
        ordered_pathway_map = OrderedDict()
        for pathway in dfa_paths:
            if pathway in pathways_reactions_map:
                ordered_pathway_map[pathway] = pathways_reactions_map[pathway]
        ordered_pathway_map = OrderedDict(reversed(list(ordered_pathway_map.items())))
        fc_plot(ax[i // 2][j], data[list(data.keys())[i]][0], ordered_pathway_map,
                condition, common_pathways, i, False)
        j += 1


    ALPHABET = ".ABCDEFGHIJKLMNOPQRSTUVWXYZ"



    ax[0][0].set_title("FASTCORE", fontsize=8)
    ax[0][1].set_title("GIMME", fontsize=8)

    plt.tight_layout()
    plt.subplots_adjust(hspace=0.40)
    plt.text(1.25, 1.15, 'NaCl Stress', fontsize=10, transform=ax[0][0].transAxes, fontweight='bold', va='top', ha='center')
    plt.text(1.25, 1.20, 'Sorbitol Stress', fontsize=10, transform=ax[1][0].transAxes, fontweight='bold', va='top', ha='center')
    plt.text(1.25, 1.20, 'H2O2 Stress', fontsize=10, transform=ax[2][0].transAxes, fontweight='bold', va='top', ha='center')
    for i in range(3):
        for j in range(2):
            ax[i][j].text(-1.5, 1.05, f"{ALPHABET[i*2 + j + 1]}", fontsize=10, transform=ax[i][j].transAxes, fontweight='bold', va='top', ha='right')
    for i in range(2):
        for j in range(2):
            ax[i][j].set_xlabel('')


    plt.savefig(f"{RESULTS_PATH}/dsalina_stress.pdf", bbox_inches='tight', format='pdf', dpi=1200)
    plt.show()
    plt.close()



if __name__ == '__main__':
    DATA_PATH = "/home/ecunha/omics-integration/data/dsalina"
    # RESULTS_PATH = "/home/ecunha/omics-integration/results/dsalina/PRJNA437866"
    RESULTS_PATH = "/home/ecunha/omics-integration/results/dsalina/PRJNA495151"
    # reaction_analysis()
    # reaction_analysis_stress()
    reaction_analysis_all()
    # gene_analysis()
