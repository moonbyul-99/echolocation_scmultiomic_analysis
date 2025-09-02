import scanpy as sc 
import matplotlib.pyplot as plt

'''
This script is used to generate the expression patter of OTOF in ZWS hippocampus additional scRNA-seq data
'''

if __name__ == '__main__':
    scdata = sc.read_h5ad('/home/rsun@ZHANGroup.local/sly_data/sly_05_require/rna_label.h5ad')
    from matplotlib.colors import LinearSegmentedColormap

    # 创建从灰色到红色的 colormap
    cmap_gray_to_red = LinearSegmentedColormap.from_list(
        'gray_to_red',
        ['lightgray', 'gray', 'darkred'],  # 可选：['gray', 'red'] 或更细腻过渡
        N=256
    )
    #sc.pl.umap(scdata, color = ['Otof'], s = 20, vmax = 1.5,cmap = 'viridis')


    ax = sc.pl.umap(scdata, color = ['Otof'], s = 20, vmax = 2,cmap = 'viridis',show =False)
    plt.savefig('OTOF.pdf', bbox_inches='tight')