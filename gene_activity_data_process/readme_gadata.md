This dir contains the source code for our gene activity data analysis. Please run the code or notebook in the order.

- 0_make_mus_gadata.ipynb is the notebook to generate the gene activity data for mouse.
- 1_make_ga is the dir contains the source code to generate the gene activity data for other 4 species
- 2_ga_data.ipynb is the notebooke to merge the 5 speices gadata and perform basic QC and analysis.
- 4_harmony_gadata.py is the script to perform harmony integration of gadata 
- 5_ga_finneanno.ipynb is the notebook to perform fine cell annotation using gadata and compare gadata annotation with scRNA-seq annotation.
- anno_info.csv is the annotation information for scrna-seq data, just the obs info of the final_anno.h5ad 

This dir is responds to the figure xx in our paper.