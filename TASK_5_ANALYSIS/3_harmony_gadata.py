import numpy as np 
import scanpy as sc 
import pandas as pd 
import matplotlib.pyplot as plt 
import scanpy.external as sce

'''
perform harmony integration among gene activity data
'''

scdata = sc.read_h5ad('/temp/sunrui/gadata_pca.h5ad')
sce.pp.harmony_integrate(scdata, 
                         key = 'species', 
                         max_iter_harmony = 50,
                         sigma =0.1, 
                         random_state = 0,
                         theta = 1)

sc.pp.neighbors(scdata, use_rep = 'X_pca_harmony')
sc.tl.umap(scdata, min_dist = 0.5)
sc.pl.umap(scdata, color = ['tissues', 'species'])

scdata.write('gadata_harmony_1.h5ad')