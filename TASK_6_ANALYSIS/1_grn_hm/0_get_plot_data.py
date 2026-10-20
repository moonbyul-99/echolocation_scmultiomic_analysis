import mudata
import numpy as np
import pandas as pd 
import matplotlib.pyplot as plt
import scanpy as sc
import anndata
from scenicplus.RSS import (regulon_specificity_scores, plot_rss)

## load scenic+ result
scplus_mdata = mudata.read("/home/rsun@ZHANGroup.local/sly_data/grn_momics/scenic_hm_outs/scenic_pipe/Snakemake/scplusmdata.h5mu")
scplus_mdata

## load scenic+ infered grn structure 
grn = scplus_mdata.uns['direct_e_regulon_metadata']
grn_e = scplus_mdata.uns['extended_e_regulon_metadata']

## load rss data

rss_tf = regulon_specificity_scores(
    scplus_mudata = scplus_mdata,
    variable = "scRNA_counts:leiden",
    modalities = ["direct_gene_based_AUC"]#, "extended_gene_based_AUC"]
)
rss_re = regulon_specificity_scores(
    scplus_mudata = scplus_mdata,
    variable = "scRNA_counts:leiden",
    modalities = ["direct_region_based_AUC"]#, "extended_gene_based_AUC"]
)
rss_idx = [str(i) for i in range(11)]
rss_tf = rss_tf.loc[rss_idx, :]
rss_re = rss_re.loc[rss_idx, :]

# modify the regulon name, uss only (+/+) or (-/-) type regulon

sel_id = []
new_col = []

for key in rss_tf.columns:
    x = key.split('_')
    if x[2] == '+/+':
        sel_id.append(True)
        new_col.append(f'{x[0]}(+)')
    elif x[2] == '-/-':
        sel_id.append(True)
        new_col.append(f'{x[0]}(-)')
    else:
        sel_id.append(False)

rss_tf = rss_tf.loc[:,sel_id]
rss_tf.columns = new_col


sel_id = []
new_col = []

for key in rss_re.columns:
    x = key.split('_')
    if x[2] == '+/+':
        sel_id.append(True)
        new_col.append(f'{x[0]}(+)')
    elif x[2] == '-/-':
        sel_id.append(True)
        new_col.append(f'{x[0]}(-)')
    else:
        sel_id.append(False)

rss_re = rss_re.loc[:,sel_id]
rss_re.columns = new_col
print(rss_re.shape)

assert (rss_re.columns == rss_tf.columns).all()
rss = rss_tf.copy() #np.sqrt(rss_re * rss_tf)
rss.to_csv('rss_new.csv')  # save rss data

## preprocessing rna data
rna = scplus_mdata['scRNA_counts']

tf_list = rss.columns.values.copy() 
for i,ele in enumerate(tf_list):
    tf_list[i] = ele.split('(')[0]
X = rna[:,tf_list]

data_list = []
for i in range(11):
    subdata = X[X.obs.leiden == str(i),:]
    data_list.append(subdata.X.mean(0))
tf_exp = pd.DataFrame(data_list, columns = rss.columns.values , index = [str(i) for i in range(11)])
tf_exp.to_csv('tf_exp.csv')  ## save tf data