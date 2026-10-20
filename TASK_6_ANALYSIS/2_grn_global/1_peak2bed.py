import numpy as np 
import pandas as pd
import muon as mu 
import scanpy as sc
import os 
import json
from multiprocessing import Pool

'''
This script is used to generate the .bed file from single cell multiomics .h5 file
'''

def get_peaks(args):
    file_name, file_path, save_dir = args
    save_path = os.path.join(save_dir, file_name)
    os.makedirs(save_path, exist_ok=True) 

    mdata = mu.read_10x_h5(file_path)
    atac = mdata['atac']

    peaks = atac.var.index.values 

    chr = []
    start = []
    end = []
    name = []

    for ele in peaks:
        a,b = ele.split(":")
        c,d = b.split("-")
        e = f'{a}_{c}_{d}'
        chr.append(a)
        start.append(c)
        end.append(d)
        name.append(e)

    with open(os.path.join(save_path, 'peaks.bed'), 'w') as f:
        for i in range(len(peaks)):
            f.write(f"{chr[i]}\t{start[i]}\t{end[i]}\t{name[i]}\n")
    print(f'{file_name} done')
    return None

def load_h5dic():
    with open('h5_dic.json', 'r') as f:
        h5_dic = json.load(f)
    
    args_list = []
    for key in h5_dic.keys():
        file_name = key 
        file_path = h5_dic[key]
        save_dir = 'peaks'
        args_list.append((file_name, file_path, save_dir))
    return args_list 

if __name__ == '__main__':
    args_list = load_h5dic()
    with Pool(processes=40) as pool:
        pool.map(get_peaks, args_list)
    print('DONE')