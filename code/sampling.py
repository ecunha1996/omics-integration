import argparse
import itertools
import os
import numpy as np
import pandas as pd
from cobra.flux_analysis.loopless import _add_cycle_free
from joblib import Parallel, delayed
from matplotlib import pyplot as plt
from scipy import stats
from scipy.stats import hypergeom, ks_2samp
from statsmodels.stats import multitest
from os import listdir
import warnings
import logging
logging.getLogger('cobra').setLevel(logging.CRITICAL)
logging.getLogger('cobra.io').setLevel(logging.CRITICAL)
warnings.filterwarnings("ignore")
from gsmmutils.model.COBRAmodel import MyModel
from cobra.io import read_sbml_model
from cobra.sampling import ACHRSampler

def split_reversible_reactions(model_to_sample):
    exchanges_demands_sinks = [reaction.id for reaction in model_to_sample.exchanges] + [reaction.id for reaction in model_to_sample.demands] + [reaction.id for reaction in model_to_sample.sinks]
    exchanges_demands_sinks = set(exchanges_demands_sinks)
    new_reactions = []
    for reaction in model_to_sample.reactions:
        if reaction not in exchanges_demands_sinks:
            if reaction.lower_bound < 0 < reaction.upper_bound:
                new_reaction = reaction.copy()
                new_reaction.id = reaction.id + "_reverse"
                new_reaction.lower_bound = 0
                new_reaction.upper_bound = -reaction.lower_bound
                for metabolite, coefficient in new_reaction.metabolites.items():
                    new_reaction.add_metabolites({metabolite: -coefficient})
                    new_reaction.add_metabolites({metabolite: -coefficient})
                new_reactions.append(new_reaction)
                reaction.lower_bound = 0
    model_to_sample.add_reactions(new_reactions)
    return model_to_sample

def load_model(filename):
    model_to_sample = MyModel(filename, "e_Biomass__cytop")
    # model_to_sample = MyModel("/home/ecunha/omics-integration/data/ngaditana/models/model_ng.xml", "e_Biomass__cytop")
    model_to_sample.objective = "e_Biomass__cytop"
    #model_to_sample.remove_reactions(find_blocked_reactions(model_to_sample), remove_orphans=True)
    for exchange in model_to_sample.exchanges:
        if exchange.lower_bound < 0 and not exchange.id.startswith("EX_C00205"):
                exchange.lower_bound = -10
        elif exchange.id.startswith("EX_C00205"):
                exchange.lower_bound = -20
    for demand in model_to_sample.demands:
        if 'photon' not in demand.id:
            demand.bounds = (0, 0)
    assert model_to_sample.optimize().status == "optimal"
    return model_to_sample#

def achr_sample_old(filename, biomass_reaction):
    try:
        model = load_model(filename)
        sol = model.optimize(objective_sense=None)
        fluxes = sol.fluxes
        opt = model.slim_optimize()
        print(opt)
        prob = model.problem
        loopless_obj_constraint = prob.Constraint(
            model.objective.expression,
            lb=opt*0.15,
            name="loopless_obj_constraint",
        )
        model.add_cons_vars([loopless_obj_constraint])
        _add_cycle_free(model, fluxes)
        solution = model.optimize(objective_sense=None)
        solution.objective_value = loopless_obj_constraint.primal
        print(solution.objective_value)
        model.remove_cons_vars(loopless_obj_constraint)
        model.reactions.e_Biomass__cytop.bounds = (opt*0.15, opt*0.6)
        sampler = ACHRSampler(model, thinning=1, seed=42)
        samples = sampler.sample(5)
        return {f"{filename.split('/')[-1].split('.xml')[0]}": samples}
    except Exception as e:
        print(e)
        print(f"Error in {filename}")

def achr_sample(filename, biomass_reaction):
    try:
        print('Sampling model: ', filename)
        model_to_sample = read_sbml_model(filename)
        for exchange in model_to_sample.exchanges:
            if exchange.lower_bound < 0 and not exchange.id.startswith("EX_C00205"):
                exchange.lower_bound = -10
            elif exchange.id.startswith("EX_C00205"):
                exchange.lower_bound = -20
        for demand in model_to_sample.demands:
            if 'photon' not in demand.id:
                demand.bounds = (0, 0)
        model_to_sample.exchanges.EX_C00244__dra.bounds = (-0.133, 1000)
        # model.exchanges.EX_C00011__dra.bounds = (-1.5, 1000)
        model_to_sample.objective = biomass_reaction
        initial_solution = model_to_sample.slim_optimize()
        print(f"Initial solution: {filename.split('/')[-1], initial_solution}")
        model_to_sample = split_reversible_reactions(model_to_sample)
        solution_after_spliting_reversible_reactions = model_to_sample.optimize().objective_value
        assert round(initial_solution, 5) == round(solution_after_spliting_reversible_reactions, 5)
        model_to_sample.reactions.get_by_id(biomass_reaction).lower_bound = initial_solution * 0.15
        solution = model_to_sample.optimize(objective_sense=None)
        print(solution.objective_value)
        sampler = ACHRSampler(model_to_sample, thinning=100, seed=42)
        samples = sampler.sample(1)
        return {f"{filename.split('/')[-1].split('.xml')[0]}": samples}
    except Exception as e:
        print(e)
        print(f"Error in {filename}")


def load_results(file_names):
    data_list = []
    for filename in file_names:
        dataframe = pd.read_csv(filename)
        if "Unnamed: 0" in dataframe.columns:
            dataframe.drop(columns=["Unnamed: 0"], inplace=True)
        # dataframe = dataframe * dataframe[biomass_map[filename.split(".")[0]]]
        data_list.append(dataframe)
    reactions = list(set([reaction for dataframe in data_list for reaction in dataframe.columns]))
    for reaction in reactions:
        for index, value in enumerate(data_list):
            if reaction not in value.columns:
                data_list[index][reaction] = 0
    results = {}
    for reaction in reactions:
        data = [value for dataframe in data_list for value in [dataframe[reaction]]]
        if not np.allclose(data[0], data[1]):
            # and not np.all_close(data[1], data[2]) and  not np.all_close(data[2], data[3]):
            stat, p_value = stats.kstest(data[0], data[1])
            print(f"Reaction: {reaction}")
            print("Statistic: ", stat)
            print("p-value: ", p_value)
            results[reaction] = [stat, p_value]
    results_df = pd.DataFrame.from_dict(results, orient="index", columns=["Statistic", "p-value"])
    results_df.to_csv("ACHR_results.csv")


def kstest(samples_control: pd.DataFrame, samples_condition: pd.DataFrame, dataset_name: str, results_dir: str):
    """
    Calculate the K-S test to detect significantly altered reactions fluxes.
    Results are saved in a csv file.

    Parameters
    ----------
    samples_control: pd.DataFrame
        The samples of the healthy tissue.
    samples_condition: pd.DataFrame
        The samples of the infected tissue.
    dataset_name: str
        The name of the dataset.
    Returns
    -------
    pd.DataFrame: The results of the K-S test for each reaction.
    :param results_dir:

    """
    print("K-S test")
    union = set(samples_condition.columns).union(set(samples_control.columns))
    samples_condition_dict, samples_control_dict = {}, {}
    for rxn in union:
        if rxn not in samples_condition.columns:
            samples_condition_dict[rxn] = pd.Series(np.zeros(10000))
        if rxn not in samples_control.columns:
            samples_control_dict[rxn] = pd.Series(np.zeros(10000))

    samples_condition = pd.concat([samples_condition, pd.DataFrame.from_dict(samples_condition_dict)], axis=1)
    samples_control = pd.concat([samples_control, pd.DataFrame.from_dict(samples_control_dict)], axis=1)

    rxns1 = set(samples_condition.columns)
    rxns2 = set(samples_control.columns)

    rxn_c = rxns1.intersection(rxns2)
    # rxn_c = rxns1.union(rxns2)
    pvals = []
    rxnid = []
    fc = []

    for rxn in rxn_c:
        data1 = samples_condition[rxn].round(decimals=4)
        data2 = samples_control[rxn].round(decimals=4)

        data1 = data1.sample(n=1000)
        data2 = data2.sample(n=1000)

        if (data1.std() != 0 and data1.mean() != 0) or (data2.std() != 0 and data2.mean() != 0):
            kstat, pval = ks_2samp(data1, data2)

            data_1_mean = data1.mean()
            data_2_mean = data2.mean()

            foldc = (data_1_mean - data_2_mean) / abs(data_1_mean + data_2_mean)

            if data_1_mean < 0 and data_2_mean < 0:
                foldc = -foldc

            pvals.append(pval)
            rxnid.append(rxn)
            fc.append(foldc)

    data_mwu = pd.DataFrame({'Reaction': rxnid, 'Pvalue': pvals})
    data_mwu = data_mwu.set_index('Reaction')

    reject, padj, _, _ = multitest.multipletests(data_mwu['Pvalue'], alpha=0.05, method='fdr_bh', is_sorted=False, returnsorted=False)

    data_mwu['Padj'] = padj
    data_mwu['Reject'] = reject
    data_mwu['FC'] = fc

    data_mwu.to_csv(f"{os.path.join(results_dir, dataset_name)}_all_results.tsv", sep="\t")

    data_sig_fc = data_mwu.loc[(abs(data_mwu['FC']) >= 0.33) & (data_mwu['Padj'] <= 0.05), :]
    data_sig_fc.to_csv(f"{os.path.join(results_dir, dataset_name)}_temp_results.tsv", sep="\t")

    rxns1 = set(samples_condition.columns)
    rxns2 = set(samples_control.columns)

    rxn_in1 = rxns1.difference(rxns2)
    rxn_in2 = rxns2.difference(rxns1)

    sigs = Parallel(n_jobs=os.cpu_count() - 2)(delayed(bootstrap_ci)(samples_condition[rx]) for rx in rxn_in1)
    act = [sigs[i][0] for i in range(len(sigs)) if sigs[i][1] == 1]
    sigs2 = Parallel(n_jobs=os.cpu_count() - 2)(delayed(bootstrap_ci)(samples_control[rx]) for rx in rxn_in2)
    rep = [sigs2[i][0] for i in range(len(sigs2)) if sigs2[i][1] == 1]

    df_abs = pd.DataFrame({'Reaction': act + rep, 'Padj_bootstrap': np.zeros(len(act + rep))})
    df_abs = df_abs.set_index('Reaction')
    # data_return = data_sig_fc + df_abs
    data_return = pd.concat([data_sig_fc, df_abs])
    file = f"{os.path.join(results_dir, dataset_name)}_results.csv"
    data_return.to_csv(file)

    return data_return.index.to_list()


def bootstrap_ci(rxn):
    """
    Calculate the confidence interval of a reaction

    Parameters
    ----------
    rxn

    Returns
    -------
    int: 1 if the reaction is significantly different from the mean, 0 otherwise.

    """
    bsci = []

    for i in range(1000):
        bt_samp = rxn.sample(1000, replace=True)
        bsci.append(bt_samp.mean())

    ci_low = np.percentile(bsci, 2.5)
    ci_high = np.percentile(bsci, 97.5)

    if ci_low > 0 or ci_high < 0:
        return rxn.name, 1
    else:
        return rxn.name, 0


def pathway_enrichment(rxnlist: list, dataset_name: str, pathways_data: pd.DataFrame, results_dir: str):
    """
    Maps significantly altered reactions to pathways using the subsystems from the HumanGEM model.
    Results are saved in csv and jpg files.

    Parameters
    ----------
    rxnlist: list
        The list of reactions to be mapped to pathways.
    dataset_name: str
        The name of the dataset.
        :param pathways_data:
    """
    try:
        print(f"Pathway enrichment for {dataset_name}...")

        listrxn_size = []
        set_size = []

        d = [g for g in rxnlist]

        for col in pathways_data.columns:
            dataframe = pd.DataFrame({'Reaction': pathways_data[col]})

            out = []

            for reac in dataframe['Reaction']:
                if reac in rxnlist:
                    out.append(reac)
                    if reac in d:
                        d.remove(reac)

            listrxn_size.append(len(out))
            set_size.append(len(pathways_data[col].dropna()))

        hyperdata = pd.DataFrame({'Pathways': pathways_data.columns, 'ListReactions': listrxn_size, 'SetSize': set_size})

        hits = hyperdata['ListReactions']
        pool = hyperdata['SetSize']

        allrxns = hyperdata['SetSize'].sum()
        targetrxns = hyperdata['ListReactions'].sum()

        pval_list = []

        for h, p in zip(hits, pool):
            rv = hypergeom(allrxns - p, p, targetrxns)

            pval = rv.pmf(h)

            pval_list.append(pval)

        hyperdata['P-value'] = pval_list

        reject, padj, _, _ = multitest.multipletests(hyperdata['P-value'], alpha=0.05, method='fdr_bh', is_sorted=False, returnsorted=False)

        hyperdata['P-value_adj'] = padj
        hyperdata['Reject'] = reject

        hyperdata_sig = hyperdata[(hyperdata['Reject']) & (hyperdata['ListReactions'] != 0)]

        hyperdata_sorted = hyperdata_sig.sort_values(by='P-value_adj', ascending=False)
        to_ignore = ['Transporters pathway', 'Biosynthesis of cofactors', 'Microbial metabolism in diverse environments']  # 'Drains pathway',
        hyperdata_sorted = hyperdata_sorted[~hyperdata_sorted['Pathways'].isin(to_ignore)]

        plt.figure(figsize=(15, 10))

        sc = plt.scatter(hyperdata_sorted['P-value_adj'], np.arange(0, len(hyperdata_sorted['Pathways'])), s=hyperdata_sorted['ListReactions'], color=(0.9, 0.3, 0.1, 0.9))

        plt.xlabel('Adjusted p-value')

        plt.yticks(np.arange(0, len(hyperdata_sorted['Pathways'])), labels=hyperdata_sorted['Pathways'])

        handles, labels = sc.legend_elements(prop="sizes", alpha=0.8)

        plt.legend(handles, labels, bbox_to_anchor=(1.6, 1.02), loc='upper right', title="Reactions")

        plt.tight_layout()

        plt.savefig(f'{os.path.join(results_dir, dataset_name)}.png', dpi=600)

        hyperdata_sorted.to_csv(f'{os.path.join(results_dir, dataset_name)}.csv', index=False)
    except Exception as e:
        print(e)


def remove_exchanges(dataframe: pd.DataFrame) -> pd.DataFrame:
    return dataframe.drop([col for col in dataframe.columns if col.startswith("EX_") or col.startswith("e_")], axis=1)


def remove_reverse_reactions(df: pd.DataFrame) -> pd.DataFrame:
    """
    Receives a dataframe with the reaction names in the columns. Removes the reactions that end with 'reverse', and sums its value to the respective reaction ID column
    Parameters
    ----------
    df

    Returns
    -------

    """
    for col in df.columns:
        if col.endswith("reverse"):
            if df[col].mean() < df[col.replace("_reverse", "")].mean():
                df[col.replace("_reverse", "")] -= df[col]
            else:
                df[col.replace("_reverse", "")] = df[col] - df[col.replace("_reverse", "")]
            df.drop(columns=[col], inplace=True)
            df[col.replace("_reverse", "")] = df[col.replace("_reverse", "")].round(6)
    return df


def read_args():
    parser = argparse.ArgumentParser(description='Process some integers.')

    # Add the arguments
    parser.add_argument('--config', type=str, required=False, help='Configuration file')
    parser.add_argument('--dataset', type=str, required=False, help='The dataset to use')
    parser.add_argument('--model', type=str, required=False, help='The model to use')
    parser.add_argument("--specific_models", type=str, required=False, help="List of filenames")
    parser.add_argument("--specific_models_dir", type=str, required=False, help="List of filenames")
    args = parser.parse_args()

    return vars(args)


if __name__ == '__main__':
    params = read_args()
    if not os.path.exists(os.path.join(os.path.dirname(params.get("model", "./")), "pathways_map.csv")):
        model = MyModel(params.get("model"), "e_Biomass__cytop")
        model.get_pathway_reactions_map()
        results_dataframe = pd.DataFrame.from_dict(data=model.pathway_reactions_map, orient='index').T
        results_dataframe.to_csv(os.path.join(os.path.dirname(params.get("model", "./")), "pathways_map.csv"), index=False)
    filenames = [os.path.join(params.get("specific_models_dir"), file) for file in listdir(params.get("specific_models_dir"))
                 if (file.endswith("3_fastcore_0.3.xml"))
                 ]
    sampling_filenames = Parallel(n_jobs=max(30, len(filenames)))(delayed(achr_sample)(filename, "e_Biomass__cytop") for filename in filenames)
    sampling_filenames = {k: v for d in sampling_filenames for k, v in d.items()}
    # sampling_filenames = {file.split("/")[-1].replace(".xml", ""):None for file in filenames}
    folder = '/'.join(filenames[0].split("/")[:-2])
    # if os.path.exists(f"{folder}/samples.h5"):
    #     os.remove(f"{folder}/samples.h5")
    for sample_name, sample in sampling_filenames.items():
        sample.to_csv(f"{folder}/{sample_name}_ACHR_samples.csv", index=False, sep="\t")
    #
    print("Loading results...")
    samples = {filename.split("_")[0] + "_" + filename.split("_")[1] : pd.read_csv(f"{folder}/{filename}_ACHR_samples.csv", sep="\t").round(6) for filename in sampling_filenames.keys()}
    # samples = {filename.split("_")[1] + "_" + filename.split("_")[-4]: None for filename in sampling_filenames}
    # for sample in samples.values():
    #     sample.drop([col for col in sample.columns if col.startswith("EX_") or col.startswith("DM_") or col.startswith("Sk_")
    #                  or col.startswith("e_")], axis=1, inplace=True)
    kstest_results = {}
    pairs = [e for e in itertools.combinations(samples, 2) if e[0].split("_")[1] == e[1].split("_")[1] and e[0].split("_")[0][-1] == e[1].split("_")[0][-1]]
    for pair in pairs:
        kstest_results['_'.join(pair)] = kstest(samples.get(pair[0]), samples.get(pair[1]), '_'.join(pair), '../results/ngaditana/PRJNA589063/dfa/')
    samples = {key: remove_reverse_reactions(sample) for key, sample in samples.items()}
    for sample_name, sample in samples.items():
        sample.to_csv(f"{folder}/{sample_name}_ACHR_samples_no_reverse.csv", index=False, sep="\t")
    dataset = pd.read_csv(os.path.join(os.path.dirname(params.get("model", "./")), "pathways_map.csv"))
    for pair, kstest_result in kstest_results.items():
        pathway_enrichment(kstest_result, pair, dataset, '../results/ngaditana/PRJNA589063/dfa/')
