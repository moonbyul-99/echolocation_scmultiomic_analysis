#module load cluster/wice/bigmem
#module load BEDTools/2.30.0-GCC-10.3.0

REGION_BED="/home/rsun@ZHANGroup.local/sly_data/sly_07_exfig/global_grn/scenic_brain_out/consensus_peak_calling/consensus_regions.bed"
GENOME_FASTA="/home/rsun@ZHANGroup.local/sly_data/genome/10x_genome/refdata-gex-mm10-2020-A/fasta/genome.fa"
CHROMSIZES="/home/rsun@ZHANGroup.local/sly_data/grn_momics/mm10.chrom.sizes"
DATABASE_PREFIX="scenic_brain_1kb_bg_with_mask"
SCRIPT_DIR="/home/rsun@ZHANGroup.local/projects_list/create_cisTarget_databases-master"

#!${SCRIPT_DIR}/create_fasta_with_padded_bg_from_bed.sh \
/home/rsun@ZHANGroup.local/projects_list/create_cisTarget_databases-master/create_fasta_with_padded_bg_from_bed.sh \
        ${GENOME_FASTA} \
        ${CHROMSIZES} \
        ${REGION_BED} \
        mm10.scenic_hm.with_1kb_bg_padding.fa \
        1000 \
        yes