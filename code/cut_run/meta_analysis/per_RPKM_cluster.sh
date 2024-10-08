# !/bin/bash
#SBATCH --job-name=computeMatrix_clusters
#SBATCH --mail-user=xxxxx
#SBATCH --mail-type=BEGIN,END
#SBATCH --account=xxxxxxx
#SBATCH --partition=largemem
#SBATCH --time=3:00:00
#SBATCH --mem=163840M

module load Bioinformatics
module load python 
source ~/miniconda3/etc/profile.d/conda.sh  # change if needed
conda activate deeptools

################## To run 4 clusters per histone mark
## H3K4me2
computeMatrix scale-regions -S H3K4me2_WT_IgG_ratio.bw -R RPKM_gene_cluster_1_WT.bed RPKM_gene_cluster_2_WT.bed RPKM_gene_cluster_3_WT.bed RPKM_gene_cluster_4_WT.bed -b 3000 -a 3000 --regionBodyLength 5000 -o RPKM_H3K4me2_WT_cluster_matrix.gz --skipZeros --outFileNameMatrix RPKM_H3K4me2_WT_cluster_scaled.tab --outFileSortedRegions RPKM_sorted_H3K4me2_WT_cluster_genes.bed
#By default, computeMatrix will handle the flipping for the negative strand, ensuring that TSS is aligned for genes on both strands in your plotProfile or plotHeatmap visualizations.
plotProfile -m RPKM_H3K4me2_WT_cluster_matrix.gz -out RPKM_H3K4me2_WT_cluster_genes_se.pdf --plotTitle "H3K4me2 distributions" --plotType=se --plotHeight 10 --plotWidth 12 --legendLocation best --plotFileFormat pdf

## H3K9me2
computeMatrix scale-regions -S H3K9me2_WT_IgG_ratio.bw -R RPKM_gene_cluster_1_WT.bed RPKM_gene_cluster_2_WT.bed RPKM_gene_cluster_3_WT.bed RPKM_gene_cluster_4_WT.bed -b 3000 -a 3000 --regionBodyLength 5000 -o RPKM_H3K9me2_WT_cluster_matrix.gz --skipZeros --outFileNameMatrix RPKM_H3K9me2_WT_cluster_scaled.tab --outFileSortedRegions RPKM_sorted_H3K9me2_WT_cluster_genes.bed
#By default, computeMatrix will handle the flipping for the negative strand, ensuring that TSS is aligned for genes on both strands in your plotProfile or plotHeatmap visualizations.
plotProfile -m RPKM_H3K9me2_WT_cluster_matrix.gz -out RPKM_H3K9me2_WT_cluster_genes_se.pdf --plotTitle "H3K9me2 distributions" --plotType=se --plotHeight 10 --plotWidth 12 --legendLocation best --plotFileFormat pdf

## H4K20me2
computeMatrix scale-regions -S H4K20me2_WT_IgG_ratio.bw -R RPKM_gene_cluster_1_WT.bed RPKM_gene_cluster_2_WT.bed RPKM_gene_cluster_3_WT.bed RPKM_gene_cluster_4_WT.bed -b 3000 -a 3000 --regionBodyLength 5000 -o RPKM_H4K20me2_WT_cluster_matrix.gz --skipZeros --outFileNameMatrix RPKM_H4K20me2_WT_cluster_scaled.tab --outFileSortedRegions RPKM_sorted_H4K20me2_WT_cluster_genes.bed
#By default, computeMatrix will handle the flipping for the negative strand, ensuring that TSS is aligned for genes on both strands in your plotProfile or plotHeatmap visualizations.
plotProfile -m RPKM_H4K20me2_WT_cluster_matrix.gz -out RPKM_H4K20me2_WT_cluster_genes_se.pdf --plotTitle "H4K20me2 distributions" --plotType=se --plotHeight 10 --plotWidth 12 --legendLocation best --plotFileFormat pdf









