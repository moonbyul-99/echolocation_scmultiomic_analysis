import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import os 
import json
from tqdm import tqdm
from scipy.sparse import isspmatrix, coo_matrix
from multiprocessing import Pool, cpu_count

'''
This script is used to convert the peak count .h5ad file to fragment file.
'''

# 全局变量（只读）
X_coo = None

def init_worker(matrix):
    global X_coo
    X_coo = matrix

def process_cell_chunk(args):
    chunk, peaks, cells = args
    results = []
    for i in chunk:
        mask = (X_coo.row == i)
        for j, count in zip(X_coo.col[mask], X_coo.data[mask]):
            chrom, start, end = peaks.iloc[j]
            results.append((chrom, int(start), int(end), cells[i], int(count)))
    return results

def adata_to_count_fragments_parallel(adata, output_dir, chr = 'chr_m', start = 'start_m', end = 'end_m', n_cpu=None):
    global X_coo

    peaks = adata.var[[chr, start, end]].copy()
    peaks[start] = peaks[start].astype(int)
    peaks[end] = peaks[end].astype(int)
    cells = adata.obs_names.tolist()

    X = adata.X #if not hasattr(adata.X, "toarray") else adata.X.toarray()
    X_coo = coo_matrix(X)

    # Step 3: 分块
    n_cpu = n_cpu or min(cpu_count(), 16)
    cell_indices = np.arange(X_coo.shape[0])
    chunks = np.array_split(cell_indices, n_cpu)

    pool_args = [(chunk, peaks, cells) for chunk in chunks]

    # Step 4: 多进程 + 进度条
    with Pool(n_cpu, initializer=init_worker, initargs=(X_coo,)) as pool:
        results = []
        for result in tqdm(pool.imap_unordered(process_cell_chunk, pool_args),
                           desc="Processing fragments",
                           total=len(pool_args)):
            results.append(result)

    # 扁平化
    fragments = [item for sublist in results for item in sublist]

    # 写入 TSV
    df = pd.DataFrame(fragments, columns=['chrom', 'start', 'end', 'cell_barcode', 'count'])
    df = df.sort_values(['chrom', 'start'])

    os.makedirs(output_dir, exist_ok=True)
    df.to_csv(os.path.join(output_dir, 'fragment.tsv'), sep='\t', index=False, header=False)
    print(f"Fragment file saved to {output_dir}")
    return df


def pipe(data_path, sample):

    #save_dir = 'fragments'
    #os.makedirs(save_dir, exist_ok=True)
    '''load data'''
    scdata = sc.read_h5ad(data_path)
    multi_data = sc.read_h5ad('/home/rsun@ZHANGroup.local/cross_embed/data/species_raw.h5ad')
    print(data_path)
    print(scdata.shape, multi_data.shape)

    '''data filter'''
    mutual_idx = np.intersect1d(multi_data.obs.index.values, scdata.obs.index.values)
    print(len(mutual_idx))

    scdata = scdata[mutual_idx,:]
    multi_data = multi_data[mutual_idx,:]

    scdata.obs = multi_data.obs.copy()

    new_idx = []
    for ele in scdata.obs.index:
        a,b,c = ele.split('-')
        x = '-'.join([b,c])
        new_idx.append(x)
    scdata.obs.index = new_idx


    '''add chr info'''
    chr = []
    start = []
    end = []

    for ele in scdata.var.index.values:
        a,_ = ele.split(':')
        b,c = _.split('-')
        chr.append(a)
        start.append(int(b))
        end.append(int(c)) 

    scdata.var['chr_m'] = chr
    scdata.var['start_m'] = start
    scdata.var['end_m'] = end

    '''save data'''
    scdata.write(data_path)

    '''atac to fragment'''
    df = adata_to_count_fragments_parallel(scdata, output_dir = f'fragments/{sample}_fragments', chr = 'chr_m', start = 'start_m', end = 'end_m', n_cpu=64)
    print(f'{sample} atac count {scdata.X.nnz}, fragment file size  {df.shape[0]}')
    print(f'{sample} DONE')

if __name__ == '__main__':
    #pipe('/home/rsun@ZHANGroup.local/sly_data/grn_momics/data/filter/cy_atac.h5ad', 'cy')
    #pipe('/home/rsun@ZHANGroup.local/sly_data/grn_momics/data/filter/jt_atac.h5ad', 'jt')
    #pipe('/home/rsun@ZHANGroup.local/sly_data/grn_momics/data/filter/qf_atac.h5ad', 'qf')
    #pipe('/home/rsun@ZHANGroup.local/sly_data/grn_momics/data/filter/t_atac.h5ad', 't')
    #pipe('/home/rsun@ZHANGroup.local/sly_data/grn_momics/data/filter/m_atac.h5ad', 'm')

    #file_dir = '/home/rsun@ZHANGroup.local/sly_data/sly_07_exfig/global_grn/modify_atac'
    #file_list = os.listdir(file_dir)

    #for file in file_list:
    #    sample = file.split('_')[0]
    #    pipe(os.path.join(file_dir,file) , sample)
    #    print('=========================================================================')

    #pipe('/home/rsun@ZHANGroup.local/sly_data/sly_07_exfig/global_grn/modify_atac/CYBQ_atac.h5ad', 'CYBQ')
    pipe('/home/rsun@ZHANGroup.local/sly_data/sly_07_exfig/global_grn/modify_atac/MHM_atac.h5ad', 'MHM')