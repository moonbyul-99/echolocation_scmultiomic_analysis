OUT_DIR="/home/rsun@ZHANGroup.local/sly_data/sly_07_exfig/global_grn/scenic_brain_out"
CBDIR="/home/rsun@ZHANGroup.local/projects_list/aertslab_motif_colleciton/v10nr_clust_public/singletons"
FASTA_FILE="/home/rsun@ZHANGroup.local/sly_data/sly_07_exfig/global_grn/mm10.scenic_hm.with_1kb_bg_padding.fa"
MOTIF_LIST="/home/rsun@ZHANGroup.local/projects_list/aertslab_motif_colleciton/v10nr_clust_public/motifs.txt"
DATABASE_PREFIX="brain_1kb_bg_with_mask"
SCRIPT_DIR="/home/rsun@ZHANGroup.local/projects_list/create_cisTarget_databases-master"

"${SCRIPT_DIR}/create_cistarget_motif_databases.py" \
    -f ${FASTA_FILE} \
    -M ${CBDIR} \
    -m ${MOTIF_LIST} \
    -o ${OUT_DIR}/${DATABASE_PREFIX} \
    --bgpadding 1000 \
    -t 20 

echo "DONE"