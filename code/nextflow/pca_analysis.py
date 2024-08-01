import plotly.express as px
import os
from gsmmutils.graphics.plot import clustermap
import matplotlib.pyplot as plt
from numpy import linspace
from sklearn.metrics import pairwise_distances
import numpy as np

def plot(components, y, xaxis_title, yaxis_title, filename, legend_title="", additional_shapes=None, annotations=None, color=None,
         width=600, height=360, **kwargs):
    if not additional_shapes:
        additional_shapes = []
    if not annotations:
        annotations = []
    if not color:
        color = y
    fig = px.scatter(x=components[:, 0], y=components[:, 1], color=color, **kwargs)
    fig.update_layout(
        xaxis_title=xaxis_title,
        yaxis_title=yaxis_title,
        margin=dict(l=20, r=20, t=20, b=20)
    )
    fig.update_layout(legend_title_text=legend_title, width=width, height=height)
    fig.update_layout(coloraxis_colorbar=dict(title='Threshold', yanchor="top", y=1.02, x=1.1, ticks="outside"))
    for shape in additional_shapes:
        fig.add_shape(
            shape
        )
    for ann in annotations:
        fig.add_annotation(
            ann
        )
    fig.show()
    fig.write_image(f"{filename}.pdf", format="pdf", scale=3)

def get_clustermap_genes_by_pathway(omics, path, omics_results=None):
    os.makedirs(path, exist_ok=True)
    omics.set_pathways_counts_by_gene()
    idx = omics.degs.index
    count_degs = omics.getmm.loc[idx]
    for pathway in omics.pathways_gene_counts.index:
        if omics.pathways_gene_counts.loc[pathway, 'Counts'] > 1:
            in_pathway = [gene_id for gene_id, value in omics.model.genes_pathways_map.items() if pathway in value]
            clustermap(count_degs.loc[count_degs.index.isin(in_pathway)], title=pathway, to_show=False, 
                       path=f"{path}/{pathway.replace(' ', '_').replace('/', '_')}.png")
            plt.close('all')

def rank_models_difference(data):
    thresholds = [e.split("_")[-1] for e in data.index]
    res = {}
    for thr in thresholds:
        data_3_fastcore_CN = data.filter(regex=f".*_{thr}", axis=0)
        # non_unique_columns = [col for col in data_3_fastcore_CN.columns if data_3_fastcore_CN[col].nunique() >1]
        # df_non_unique = data_3_fastcore_CN[non_unique_columns]
        hamming_distances = pairwise_distances(data_3_fastcore_CN.astype(bool).values, metric='hamming')
        mean_hamming_distance = np.mean(hamming_distances)
        res[thr] = mean_hamming_distance
    return res