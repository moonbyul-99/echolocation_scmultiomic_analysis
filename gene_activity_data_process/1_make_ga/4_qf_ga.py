import os
import numpy as np 
import pandas as pd 
import scanpy as sc 
import anndata as ad
import snapatac2 as snap 
import json 
import logging

# 配置日志记录器
logging.basicConfig(filename='qf_process.log', level=logging.INFO, 
                    format='%(asctime)s - %(levelname)s - %(message)s')

path_dic = {
    'bq':'/home/rsun@ZHANGroup.local/sly_data/data/QF/QFBQ/atac_fragments.tsv.gz',
    'ct':'/home/rsun@ZHANGroup.local/sly_data/data/QF/QFCT/atac_fragments.tsv.gz',
    'hm':'/home/rsun@ZHANGroup.local/sly_data/data/QF/QFHM/outs/atac_fragments.tsv.gz',
    'ng':'/home/rsun@ZHANGroup.local/sly_data/data/QF/QFNG/atac_fragments.tsv.gz',
    'qn':'/home/rsun@ZHANGroup.local/sly_data/data/QF/QFQN/atac_fragments.tsv.gz',
    'pc':'/home/rsun@ZHANGroup.local/sly_data/data/QF/QFPC/atac_fragments.tsv.gz',
    'xq':'/home/rsun@ZHANGroup.local/sly_data/data/QF/QFXQ/atac_fragments.tsv.gz',
    'yb':'/home/rsun@ZHANGroup.local/sly_data/data/QF/QFYB/atac_fragments.tsv.gz'
    }

SPECIES = 'qf'

if __name__ == '__main__': 

    save_dir = f'{SPECIES}_atac_data'
    os.makedirs(save_dir, exist_ok=True)

    chrom_dic_path = '/home/rsun@ZHANGroup.local/sly_data/genome/10x_genome/csp/star/chrNameLength.txt'

    chrom_sizes = {}
    with open(chrom_dic_path, 'r') as f:
        for line in f:
            chr, length = line.strip().split('\t')
            chrom_sizes[chr] = int(length)

    ## load data
    for key in path_dic:
        if key in ['bq','ct']:
            continue
        logging.info(f'begin processing {key} data')
        frag_path = path_dic[key]
        save_path = os.path.join(save_dir, f'{SPECIES}{key}_raw.h5ad')
        scdata = snap.pp.import_data(frag_path,
                                        chrom_sizes = chrom_sizes, 
                                        sorted_by_barcode = False,
                                        file = save_path)
        
        logging.info(scdata.close(),scdata.is_closed())

        # load data into memory, avoid modifying raw data
        scdata = snap.read(save_path, backed = None)

        #calculate and plot the size distribution of fragments in this dataset
        snap.pl.frag_size_distr(scdata, interactive=False)
        gtf_file = '/home/rsun@ZHANGroup.local/sly_data/genome/10x_genome/csp/genes/genes.gtf'
        snap.metrics.tsse(scdata, gene_anno = gtf_file)
        snap.pl.tsse(scdata, interactive=False)


        logging.info(f'raw scdata shape {scdata.shape}')
        snap.pp.filter_cells(scdata, min_counts=5000, min_tsse=5, max_counts=100000)
        logging.info(f'filtered scdata shape {scdata.shape}') 

        # get gadata 

        new_gtf = '/home/rsun@ZHANGroup.local/sly_data/notebook/csp_new.gtf'
        gadata = snap.pp.make_gene_matrix(scdata, gene_anno=new_gtf)

        ## get gene symbol 

        with open('/home/rsun@ZHANGroup.local/sly_data/genome/csp_gene_map.json', 'r') as f:
            gene_map = json.load(f)

        gadata.var.loc[:,'gene_symbol'] = gadata.var.index.values 

        for i, ele in enumerate(gadata.var.loc[:,'gene_symbol'].values):
            if ele in gene_map:
                gadata.var.iloc[i,:]['gene_symbol'] = gene_map[ele].capitalize()

        ## save gadata 
        gadata.write(os.path.join(save_dir, f'{SPECIES}{key}_gadata.h5ad'))