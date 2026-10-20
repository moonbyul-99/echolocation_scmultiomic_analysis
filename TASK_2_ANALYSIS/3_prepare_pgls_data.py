import numpy as np 
import pandas as pd 
import anndata as ad 
import scanpy as sc 
# import ete3
import matplotlib.pyplot as plt 
import seaborn as sns

if __name__ == '__main__':
    scdata = sc.read_h5ad('../../sly_data/data/final_anno.h5ad')
    save_dir = 'pgls_data'
    import os 

    for tissue in scdata.obs.tissues.unique():
        res = {}
        for species in scdata.obs.species.unique():
            idx = (scdata.obs.tissues == tissue) & (scdata.obs.species == species)
            sub = scdata[idx].obsm['X_pca_harmony'][:, :] 
            n = len(sub)
            sample = sub[np.random.choice(n, size=min(n, 2000),replace=False)]
            res[species] = sample.mean(axis=0)
        df = pd.DataFrame(res).T 
        df.columns = [f'PC{i}' for i in range(1,51)]
        df.to_csv(os.path.join(save_dir, f'expr_harmony_{tissue}.csv'))

    scdata = scdata[scdata.obs.tissues == 'HM',:]
    for celltype in scdata.obs.celltype.unique():
        res = {}
        for species in scdata.obs.species.unique():
            if celltype == 'inhibitory neuron':
                idx = np.load('pgls_data/sel_idx_inhibitory_neuron_ga.npz')[species]
            elif celltype == 'excitatory neuron':
                idx = np.load('pgls_data/sel_idx_excitatory_neuron_ga.npz')[species]
            else:
                idx = (scdata.obs.celltype == celltype) & (scdata.obs.species == species)
            sub = scdata[idx].obsm['X_pca_harmony'][:, :] 
            n = len(sub)
            sample = sub[np.random.choice(n, size=min(n, 1000),replace=False)]
            res[species] = sample.mean(axis=0)
        df = pd.DataFrame(res).T 
        df.columns = [f'PC{i}' for i in range(1,1+df.shape[1])]
        df.to_csv(os.path.join(save_dir, f'expr_cell_{celltype}.csv'))