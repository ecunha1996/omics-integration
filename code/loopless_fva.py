import pandas as pd
from gsmmutils.model.COBRAmodel import MyModel
from cobra.flux_analysis import flux_variability_analysis, loopless_solution


def load_model():
    model_to_sample = MyModel("/home/ecunha/omics-integration/results/ngaditana/PRJNA589063/integration/models/N5_fastcore_1.2.xml", "e_Biomass__cytop")
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





if __name__ == "__main__":
    model = load_model()
    # fva_sol = model.loopless_fva(fraction_of_optimum=0.15, n_jobs = 50)
    # fva_sol.to_csv("/home/ecunha/omics-integration/results/ngaditana/PRJNA589063/integration/loopless_fva.tsv", sep="\t")
    fva_sol = pd.read_csv(r"/home/ecunha/omics-integration/results/ngaditana/PRJNA589063/integration/loopless_fva.tsv", sep="\t", index_col=0)
    print(loopless_solution(model))
    reactions_bounds = fva_sol.to_dict(orient='index')
    for reaction in model.reactions:
        print(reaction.id)
        if reactions_bounds[reaction.id]["maximum"]<reactions_bounds[reaction.id]["minimum"]:
            pass
        else:
            reaction.bounds = (reactions_bounds[reaction.id]["minimum"], reactions_bounds[reaction.id]["maximum"])
            sol =  model.optimize()
            if sol.status != "optimal":
                print(sol.objective_value)
                print(reaction.bounds)
                print("-"*50)
    print(model.optimize())
