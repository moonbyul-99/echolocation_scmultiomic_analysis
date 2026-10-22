## Analysis Pipeline for Convergent Evolution of a Novel Hippocampal Cell Type Underlying Mammalian Echolocation

## TASK_1_ANALYSIS (Basic scRNA-seq analysis)

This dir contains the code for reproducing the following figures:

$\color{red}{\text{Fig 1 B-D + Fig 2 AB}}$
1. Prepare data, run 1_get_h5ad.py and 2_total_data.py.
2. Perform the RNA-seq based cell annotation, including batch effect using harmony.
   Refer the 3_anno_alldata.ipynb and the 4_anno_alldata.ipynb, generate $\color{red}{\text{Figure 1B}}$
3. Plot the marker dotplot and DEG heatmap using 05_marker.ipynb and 06_conserved_marker.ipynb, generate $\color{red}{\text{Figure 1 C,D}}$
4. Using the annotation information, readers could get the data for the cell porpotion
   plot and generate $\color{red}{\text{Figure 2 A,B}}$, please refer the 07_cell_ratio.ipynb for $\color{red}{\text{Figure 2B}}$. 

## TASK_2_ANALYSIS (scRNA-seq based divergence score analysis)

This dir contains the code for reproducing the following figures: 
$\color{red}{\text{Fig 2 CD + Supp Fig 4}}$.

1. For the tissue level rna-divergence, using 1_tissue_diff.ipynb, generate $\color{red}{\text{Fig 2C}}$
2. For the celltype level rna-divergence, using 2_celltype_diff.ipynb, generate $\color{red}{\text{Fig 2D}}$
3. For the D-PGLS analysis, extracting the 
   required data from single cell data and 
   use D-PGLS package to generate $\color{red}{\text{Supp Fig 4}}$. 

## TASK_3_ANALYSIS (HM scRNA-seq analysis)

This dir contains the code for reproducing the following figures:
$\color{red}{\text{Figure 3 A-C, F}}$ + $\color{red}{\text{Supp Fig 5 A, B, C}}$ + $\color{red}{\text{Supp Fig 6}}$ + $\color{red}{\text{Supp Fig 9}}$

For creating gene-activity data from scATAC-seq data, ref the make_ga/ dir
1. HM scrna-seq fine analysis, ref 00_hmdata_full.ipynb + 01_cluster_echo_ex.ipynb, Generate $\color{red}{\text{Figure 3 A and Supp Fig 5 A}}$. The annotation result is listed in HM_cell_anno.csv.
2. Cluster resolution analysis, ref 02_cluster_check.ipynb, get the cluster resolution test and generate Supp $\color{red}{\text{Fig 5 B}}$, 
   and with the defined resolution, get the cell proportion plot, $\color{red}{\text{Figure 3C + Supp Fig 5 C}}$
3. Using the GBDT based marker gene set identification, generate the $\color{red}{\text{Figure 3B}}$ using 03_gbdt_marker.ipynb and $\color{red}{\text{Supp Fig 6}}$ using 04_modify_umap.ipynb
4. Using the 05_cluster_tree.ipynb for the $\color{red}{\text{Figure 3F}}$ 
5. For additional species analysis, ref 07_cluster_recheck.ipynb, generate the $\color{red}{\text{Supp Fig 9}}$.
We listed the data links of supplementary species  
in supplementary meterials. Some data require 
re-analysis before 07_cluster_recheck.ipynb and 
we will release these preprocessing scripts soon.

## TASK_4_ANALYSIS (Additional ZWS RNA data analysis)
This dir contains the code for reproducing the following figures:
$\color{red}{\text{Figure 4 B-D}}$ + $\color{red}{\text{Supp Fig 5 D, E}}$

1. 1_analysis.ipynb is the basic scRNA-seq analysis for the new ZWS data. 
2. 2_neuron_label.ipynb is the label transfer notebook for the ZWS rna data and generate $\color{red}{\text{figure 4B}}$. 
3. 3_rna_diff.ipynb generate the cluster differences result of $\color{red}{\text{Figure 4D}}$.
4. Analysis above directly generate the data required for $\color{red}{\text{Supp Fig 5 D, E and Figure 4C}}$.


## TASK_5_ANALYSIS (Baisc ATAC-seq analysis + QC plot)
This dir contains the code for reproducing the following figures:
$\color{red}{\text{Supplementary Figure 3A-C}}$.

1. Get gene activity data, using 2_gadata.ipynb
2. Using Harmonoy batch correction, using 3_harmony_gadata.py 
3. Perform cell type annotation and check the annotation with scRNA-seq result, using 4_gadata_annotation.ipynb + 5_ga_fineanno.ipynb, 
   generate the result of $\color{red}{\text{Supplementary Figure 3A-C}}$. $\color{red}{\text{Supp Fig 3 B}}$ is similiar with $\color{red}{\text{Figure 1 C}}$, just change the data from rna-seq to gene activity,
   so check the plot code of $\color{red}{\text{Figure 1C}}$.

## TASK_6_ANALYSIS (scenic+ GRN analysis)
This dir contains the code for reproducing the following figures:
$\color{red}{\text{Figure 1E, Figure 2E,F, Figure 3G-I,}}$ $\color{red}{\text{Supplementary Figure 3D and Supplementary Figure 1, 2}}$

1. For the brain regulation inference, please refer the 2_grn_global dir. Please read the readme_grn_global.md file for the cross species grn processing. 
   After the grn inference, Get the plot of $\color{red}{\text{Supplementary Figure 3D and Figure 1E}}$ using 11_get_species_reg_data.ipynb. Noted that the Figure 1E is a subset 
   of Supplementary Figure 3D. 
2. For the regulation difference quantification across brain regions in $\color{red}{\text{Figure 2E}}$, we generate the regulon differences file in 12_reg_outlier.ipynb. The 
   generated abs_diff_df_z.csv is the input data for the outlier analysis in Figure 2E. The Figure 2E data is calculated using outlier detection method, 
   one class SVM. For the regulon difference plot across species in HM, also refer 12_reg_outlier.ipynb.
3. For the HM specific regulation inference, please refer the 1_grn_hm dir. Noted that the HM specific GRN inference is similar with the whole brain GRN 
   inference, so we temporarily skipped the whole pipeline code and just provide the plot code. The complete code we will provide later.
   -  $\color{red}{\text{Figure 3G}}$ code, first run the 0_get_plot_data.py and then run the 1_regulon_dot_heatmap.ipynb 
   -  $\color{red}{\text{Figure 3H}}$ code, see 1_grn_hm/Figure_3H, first run the get_data.ipynb to get the required data, then using CytoScape to plot. For the GSEA analysis, 
   using the target gene in the selected eRegulons.
   -  $\color{red}{\text{Figure 3I}}$ code, see 1_grn_hm/Figure_3I, run the get_data.ipynb to get the required data, then using CytoScape to plot. For the GSEA analysis, using the 
   target genes in the selected eRegulons.
4. For $\color{red}{\text{Supplementary Figure 1, 2}}$, which showing the joint QC result in multi-omics analysis, plot code is listed in qc_plot.ipynb.
   The source of data required is listed in this notebook.
