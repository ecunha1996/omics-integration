import pandas as pd
from os.path import join
from scipy.stats import zscore


def load_samples(DATA_PATH, methods = {'fastcore': 0.6, 'gimme': 0.4}, samples = []):
    data = {}
    for method, threshold in methods.items():
        data_list = {}
        tmp = []
        for sample in samples:
            data_list[sample] = pd.read_csv(join(DATA_PATH, f"samples/{sample}_{method}_{threshold}_ACHR_samples.csv"), sep="\t")
        for name, ds in data_list.items():
            tmp.append(ds.mean())
        data_df = pd.DataFrame(tmp, index=samples).T
        data_df = data_df.fillna(0)
        row_std = data_df.std(axis=1)
        # data_df = data_df[row_std != 0]
        data_df = data_df.add_suffix(f"_{method}")
        data[method] = data_df
    total_data = pd.concat([e.apply(lambda x: zscore(x), axis=1) for e in data.values()],axis=1)
    total_data = total_data.fillna(0)
    row_std = total_data.std(axis=1)
    data['all'] = total_data[row_std != 0]
    return data