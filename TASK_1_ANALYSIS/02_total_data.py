### get all 40 h5ad into one h5ad file 
''' merge all 40 .h5ad files into one .h5ad file.'''
import numpy as np 
import pandas as pd
import matplotlib.pyplot as plt 
import seaborn as sns
import os
import scanpy as sc
import logging
import argparse
from datetime import datetime 
import json
import anndata as ad
from scipy.sparse import csr_matrix
import scanpy.external as sce
from sklearn.metrics.pairwise import cosine_similarity

if __name__ == '__main__':
    cydata = sc.read_h5ad('/work/sunrui/sly_data/ori_data/CYBQ.h5ad')
    jtdata = sc.read_h5ad('/work/sunrui/sly_data/ori_data/JTBQ.h5ad')
    qfdata = sc.read_h5ad('/work/sunrui/sly_data/ori_data/QFBQ.h5ad')
    mdata = sc.read_h5ad('/work/sunrui/sly_data/ori_data/MBQ.h5ad')
    tdata = sc.read_h5ad('/work/sunrui/sly_data/ori_data/TBQ.h5ad')
    
    mutual_gene = np.intersect1d(cydata.var.index.values, jtdata.var.index.values)
    mutual_gene = np.intersect1d(qfdata.var.index.values, mutual_gene)
    mutual_gene = np.intersect1d(mdata.var.index.values, mutual_gene)
    mutual_gene = np.intersect1d(tdata.var.index.values, mutual_gene)
    
    print('====='*20)
    print(mutual_gene.shape)
    print('====='*20)
    
    species = ['CY','JT','QF','M','T']
    tissues = ['BQ', 'CT', 'HM', 'NG', 'PC', 'QN', 'XQ', 'YB']
    data_list = []
    for s in species:
        for t in tissues:
            scdata = sc.read_h5ad(f'work/sunrui/sly_data/ori_data/{s+t}.h5ad')
            scdata = scdata[:,mutual_gene]
            data_list.append(scdata)
    final_data = ad.concat(data_list)
    print(final_data.shape)
    final_data.write_h5ad(f'/work/sunrui/sly_data/all_data.h5ad')