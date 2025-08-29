### get all 40 sample scRNA-seq h5ad file 
'''
read 40 scmultiome data(5 species, 8 tissues), rna-seq only
for each sample, for its gene_ids, if has orthologe mus gene, use ortgologe gene id, else remain the species specific gene id.
save altsl 40 .h5ad files
'''
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


##############################################################

def read_data(data_path,species,tissues):
    scdata = sc.read_10x_mtx(data_path, var_names = 'gene_symbols',cache= True)  
    scdata.obs.loc[:,'species'] = species 
    scdata.obs.loc[:,'tissues'] = tissues
    scdata.obs.loc[:,'samples'] = species + tissues
    return scdata

def drop_duplicates(adata):
    # 有些数据存在重复的基因注释，将这些重复的基因求和作为新的基因
    # some genes are duplicated, merge these duplicates as one gene
    column_name = adata.var.index.values
    for i, value in enumerate(column_name):
        column_name[i] = value.capitalize()
    tmp_adata = pd.DataFrame(adata.X.todense(), columns=column_name,index= adata.obs.index.values.reshape(-1))
    tmp_adata = tmp_adata.groupby(tmp_adata.columns,axis=1).sum()
    
    
    tmp_X = csr_matrix(tmp_adata.values)
    
    new_data = ad.AnnData(X = tmp_X)
    new_data.var.index = tmp_adata.columns
    new_data.var['gene_ids'] = new_data.var.index.values.copy()
    new_data.var['feature_types'] = 'Gene Expression'
    new_data.obs = adata.obs
    return new_data

def process_query_genome(species, query):
    # for changyi bat genome
    if species == 'CY':
        with open("/work/sunrui/sly_data/genome_transform/mfu_gene_map.json",'r') as f:
            cy_gene_map = json.load(f)
        f.close()

        query.var.loc[:,'gene_ids'] = query.var.loc[:,'gene_ids'].replace(cy_gene_map)
        query.var.index = query.var.loc[:,'gene_ids'].values.copy()
        query = drop_duplicates(query)
        
    elif species == 'JT':
        # for jutou bat genome
        with open("/work/sunrui/sly_data/genome_transform/rsi_gene_map.json",'r') as f:
            jt_gene_map = json.load(f)
        f.close()
        
        query.var.loc[:,'gene_ids'] = query.var.loc[:,'gene_ids'].replace(jt_gene_map)
        query.var.index = query.var.loc[:,'gene_ids'].values.copy()
        query = drop_duplicates(query)
        
    elif species == 'QF':
        # for quanfu bat genome 
        qf_map=pd.read_excel("/work/sunrui/sly_data/genome_transform/csp_gene.xlsx",sheet_name='CSP')
        qf_dic = {}
        for i in range(qf_map.shape[0]):
            qf_dic[qf_map.iloc[i,1]] = qf_map.iloc[i,0]
        
        query.var.loc[:,'gene_ids'] = query.var.loc[:,'gene_ids'].replace(qf_dic)
        query.var.index = query.var.loc[:,'gene_ids'].values.copy()
        query = drop_duplicates(query)
    
    elif species == 'M':
        pass
    elif species =='T':
        pass
    else:
        print('Wrong Species')
    return query

#################################################################################################################
if __name__ == '__main__':
    data_path_dict = {
        'CYBQ':'/work/sunrui/rna_seq/clean_rna/CY/result/CYBQ/outs/filtered_feature_bc_matrix',
        'CYCT':'/work/sunrui/rna_seq/clean_rna/CY/CYCT/outs/filtered_feature_bc_matrix',
        'CYHM':'/work/sunrui/rna_seq/clean_rna/CY/CYHM/outs/filtered_feature_bc_matrix',
        'CYNG':'/work/sunrui/rna_seq/clean_rna/CY/result/CYNG/outs/filtered_feature_bc_matrix',
        'CYPC':'/work/sunrui/rna_seq/clean_rna/CY/result/CYPC/outs/filtered_feature_bc_matrix',
        'CYQN':'/work/sunrui/rna_seq/clean_rna/CY/result/CYQN/outs/filtered_feature_bc_matrix',
        'CYXQ':'/work/sunrui/rna_seq/clean_rna/CY/result/CYXQ/outs/filtered_feature_bc_matrix',
        'CYYB':'/work/sunrui/rna_seq/clean_rna/CY/result/CYYB/outs/filtered_feature_bc_matrix',

        'JTBQ':'/work/sunrui/rna_seq/clean_rna/JT/result/JTBQ/outs/filtered_feature_bc_matrix',
        'JTCT':'/work/sunrui/rna_seq/clean_rna/JT/result/JTCT/outs/filtered_feature_bc_matrix',
        'JTHM':'/work/sunrui/rna_seq/clean_rna/JT/result/JTHM/outs/filtered_feature_bc_matrix',
        'JTNG':'/work/sunrui/rna_seq/clean_rna/JT/result/JTNG/outs/filtered_feature_bc_matrix',
        'JTPC':'/work/sunrui/rna_seq/clean_rna/JT/result/JTPC/outs/filtered_feature_bc_matrix',
        'JTQN':'/work/sunrui/rna_seq/clean_rna/JT/result/JTQN/outs/filtered_feature_bc_matrix',
        'JTXQ':'/work/sunrui/rna_seq/clean_rna/JT/result/JTXQ/outs/filtered_feature_bc_matrix',
        'JTYB':'/work/sunrui/rna_seq/clean_rna/JT/result/JTYB/outs/filtered_feature_bc_matrix',

        'QFBQ':'/work/sunrui/rna_seq/clean_rna/QF/result/QFBQ/outs/filtered_feature_bc_matrix',
        'QFCT':'/work/sunrui/rna_seq/clean_rna/QF/result/QFCT/outs/filtered_feature_bc_matrix',
        'QFHM':'/work/sunrui/rna_seq/clean_rna/QF/result/QFHM/outs/filtered_feature_bc_matrix',
        'QFNG':'/work/sunrui/rna_seq/clean_rna/QF/result/QFNG/outs/filtered_feature_bc_matrix',
        'QFPC':'/work/sunrui/rna_seq/clean_rna/QF/result/QFPC/outs/filtered_feature_bc_matrix',
        'QFQN':'/work/sunrui/rna_seq/clean_rna/QF/result/QFQN/outs/filtered_feature_bc_matrix',
        'QFXQ':'/work/sunrui/rna_seq/clean_rna/QF/result/QFXQ/outs/filtered_feature_bc_matrix',
        'QFYB':'/work/sunrui/rna_seq/clean_rna/QF/result/QFYB/outs/filtered_feature_bc_matrix',

        'MBQ':'/temp/sunrui/swap_backup/dataset/mus_matrix/MBQ/outs/filtered_feature_bc_matrix',
        'MCT':'/work/sunrui/rna_seq/clean_rna/M/MCT/outs/filtered_feature_bc_matrix',
        'MHM':'/temp/sunrui/swap_backup/dataset/mus_matrix/MHM/outs/filtered_feature_bc_matrix',
        'MNG':'/temp/sunrui/swap_backup/dataset/mus_matrix/MNG/outs/filtered_feature_bc_matrix',
        'MPC':'/temp/sunrui/swap_backup/dataset/mus_matrix/MPC/outs/filtered_feature_bc_matrix',
        'MQN':'/temp/sunrui/swap_backup/dataset/mus_matrix/MQN/outs/filtered_feature_bc_matrix',
        'MXQ':'/temp/sunrui/swap_backup/dataset/mus_matrix/MXQ/outs/filtered_feature_bc_matrix',
        'MYB':'/temp/sunrui/swap_backup/dataset/mus_matrix/MYB/outs/filtered_feature_bc_matrix',

        'TBQ':'/temp/sunrui/swap_backup/dataset/tda_matrix/TBQ/outs/filtered_feature_bc_matrix',
        'TCT':'/temp/sunrui/swap_backup/dataset/tda_matrix/TCT/outs/filtered_feature_bc_matrix',
        'THM':'/temp/sunrui/swap_backup/dataset/tda_matrix/THM/outs/filtered_feature_bc_matrix',
        'TNG':'/temp/sunrui/swap_backup/dataset/tda_matrix/TNG/outs/filtered_feature_bc_matrix',
        'TPC':'/temp/sunrui/swap_backup/dataset/tda_matrix/TPC/outs/filtered_feature_bc_matrix',
        'TQN':'/temp/sunrui/swap_backup/dataset/tda_matrix/TQN/outs/filtered_feature_bc_matrix',
        'TXQ':'/temp/sunrui/swap_backup/dataset/tda_matrix/TXQ/outs/filtered_feature_bc_matrix',
        'TYB':'/temp/sunrui/swap_backup/dataset/tda_matrix/TYB/outs/filtered_feature_bc_matrix'
        }
    
    species = ['CY', 'JT', 'QF', 'M', 'T']
    tissues = ['BQ', 'CT', 'HM', 'NG', 'PC', 'QN', 'XQ', 'YB']
    
    save_dir = '/work/sunrui/sly_data/ori_data/'
    if not os.path.exists(save_dir):
        os.makedirs(save_dir)
    
    
    for s in species:
        for t in tissues:
            data_path = data_path_dict[s+t]
            
            scdata = read_data(data_path,s,t)
            scdata = process_query_genome(s, scdata)
            scdata = drop_duplicates(scdata)
            scdata.write_h5ad(os.path.join(save_dir, f'{s+t}.h5ad'))
            print(f'{s+t} over, data size is {scdata.shape}')
    print('over' + '******'*10)
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    