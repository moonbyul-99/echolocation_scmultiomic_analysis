This dir contains the source code for the gene regulation network analysis in the whole brain. 

The basic pipeline is as follows:

1. get the .bed file of the scATAC-seq data. (1_peak2bed.py)
2. transfer the species specific peak to the mm10 peak. (2_peak_trans.py)
3. based on the peak transfer result, extract the transfer info into a mapping file. For each species, we construct a dictionary, records the corresponding relationship between species peaks and mm10 peaks. (3_global_peak_map.ipynb)
4. modify the rna file to the cistopic required format. (4_create_multi_rna.ipynb)
5. generate the new peak count matrix for 5 species, now the var info is their mm10 peak. (5_atac_modify.ipynb)
6. transform the peak count matrix to the fragment.tsv file require for cistopic. (6_1_*, 6_2_*, 6_3_*) 

    <font color = 'red'> 
    Noted that due to the specific format of the cell barcode required by cistopic, 6_2_* 6_3_* are used to modify some error case caused by 6_1_*. 

    A recommended cell barcode format is like {10x cell barcoder}_{sample name}. Cistopic document does not provide a clear illustration of the barcode format, however, illegal format will cause error when running cistopic.
    </font> 

7. compress the tsv file to .tsv.gz. (7_1_*, 7_2_*)
8. run the cistopic pipeline (8_scenci_brain.ipynb) 
   
   <font color = 'red'>
    Noted that running this notebook will cost a lot of time and requires a large amount of memory. Cistopic is a powerful tool, but some implementations are not very efficient. For example, it uses dense matrix instead of sparse matrix in computing, which greatly decreases the speed.
   </font>

9. cistarget pipeline (9_1_*, 9_2_*)
10. After the pipeline above, we are able to run scienplut pipeline, the Snakemake config file is in the cofig dir. 
11. grn analysis and visualization (10_* , 11_*) 

This dir is used to generate the figure xxxx in our paper.