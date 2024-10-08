# !/bin/bash
#SBATCH --job-name=computeMatrix_WT_KO
#SBATCH --mail-user=xxxx
#SBATCH --mail-type=BEGIN,END
#SBATCH --account=xxxxx
#SBATCH --partition=largemem
#SBATCH --time=18:00:00
#SBATCH --mem=163840M

module Bioinformatics
module load python 
source ~/miniconda3/etc/profile.d/conda.sh  # change if needed
conda activate deeptools


#### All genome
## H3K4me2
computeMatrix reference-point -S H3K4me2_WT_IgG_ratio.bw H3K4me2_KO_IgG_ratio.bw -R 2024_09_25_all_genes.bed -b 30000 -a 30000 --referencePoint TSS -o all_2024_09_25_H3K4me2_matrix.gz --skipZeros --outFileNameMatrix all_2024_09_25_H3K4me2_scaled.tab --outFileSortedRegions all_2024_09_25_H3K4me2_sorted.bed
#By default, computeMatrix will handle the flipping for the negative strand, ensuring that TSS is aligned for genes on both strands in your plotProfile or plotHeatmap visualizations.
plotProfile -m all_2024_09_25_H3K4me2_matrix.gz -out all_2024_09_25_H3K4me2.pdf --plotTitle "H3K4me2 distributions" --plotType=se --yMin=0.99 --yMax=1.2 --perGroup --plotHeight 10 --plotWidth 12 --legendLocation best --plotFileFormat pdf
## H3K9me2
computeMatrix reference-point -S H3K9me2_WT_IgG_ratio.bw H3K9me2_KO_IgG_ratio.bw -R 2024_09_25_all_genes.bed  -b 30000 -a 30000 --referencePoint TES -o all_2024_09_25_H3K9me2_matrix.gz --skipZeros --outFileNameMatrix all_2024_09_25_H3K9me2_scaled.tab --outFileSortedRegions all_2024_09_25_H3K9me2_sorted.bed
#By default, computeMatrix will handle the flipping for the negative strand, ensuring that TSS is aligned for genes on both strands in your plotProfile or plotHeatmap visualizations.
plotProfile -m all_2024_09_25_H3K9me2_matrix.gz -out all_2024_09_25_H3K9me2.pdf --plotTitle "H3K9me2 distributions" --plotType=se --yMin=0.99 --yMax=1.2 --perGroup --plotHeight 10 --plotWidth 12 --legendLocation best --plotFileFormat pdf
## H4K20me2
computeMatrix reference-point -S H4K20me2_WT_IgG_ratio.bw H4K20me2_KO_IgG_ratio.bw -R 2024_09_25_all_genes.bed -b 30000 -a 30000 --referencePoint TES -o all_2024_09_25_H4K20me2_matrix.gz --skipZeros --outFileNameMatrix all_2024_09_25_H3K9me2_scaled.tab --outFileSortedRegions all_2024_09_25_H4K20me2_sorted.bed
#By default, computeMatrix will handle the flipping for the negative strand, ensuring that TSS is aligned for genes on both strands in your plotProfile or plotHeatmap visualizations.
plotProfile -m all_2024_09_25_H4K20me2_matrix.gz -out all_2024_09_25_H4K20me2.pdf --plotTitle "H4K20me2 distributions" --plotType=se --yMin=0.90 --yMax=1.2 --perGroup --plotHeight 10 --plotWidth 12 --legendLocation best --plotFileFormat pdf


##### UP DEG
## H3K4me2
computeMatrix reference-point -S H3K4me2_WT_IgG_ratio.bw H3K4me2_KO_IgG_ratio.bw -R 2024_06_19_UP_DEG.bed -b 30000 -a 30000 --referencePoint TSS -o 24_06_UP_H3K4me2_matrix.gz --skipZeros --outFileNameMatrix 24_06_UP_H3K4me2_scaled.tab --outFileSortedRegions 24_06_UP_H3K4me2_sorted.bed
#By default, computeMatrix will handle the flipping for the negative strand, ensuring that TSS is aligned for genes on both strands in your plotProfile or plotHeatmap visualizations.
plotProfile -m 24_06_UP_H3K4me2_matrix.gz -out 24_06_UP_H3K4me2.pdf --plotTitle "H3K4me2 distributions" --plotType=se --yMin=0.8 --yMax=1.9 --perGroup --plotHeight 10 --plotWidth 12 --legendLocation best --plotFileFormat pdf
## H3K9me2
computeMatrix reference-point -S H3K9me2_WT_IgG_ratio.bw H3K9me2_KO_IgG_ratio.bw -R 2024_06_19_UP_DEG.bed  -b 30000 -a 30000 --referencePoint TES -o 24_06_UP_H3K9me2_matrix.gz --skipZeros --outFileNameMatrix 24_06_UP_H3K9me2_scaled.tab --outFileSortedRegions 24_06_UP_H3K9me2_sorted.bed
#By default, computeMatrix will handle the flipping for the negative strand, ensuring that TSS is aligned for genes on both strands in your plotProfile or plotHeatmap visualizations.
plotProfile -m 24_06_UP_H3K9me2_matrix.gz -out 24_06_UP_H3K9me2.pdf --plotTitle "H3K9me2 distributions" --plotType=se --yMin=0.8 --yMax=1.8 --perGroup --plotHeight 10 --plotWidth 12 --legendLocation best --plotFileFormat pdf
## H4K20me2
computeMatrix reference-point -S H4K20me2_WT_IgG_ratio.bw H4K20me2_KO_IgG_ratio.bw -R 2024_06_19_UP_DEG.bed -b 30000 -a 30000 --referencePoint TES -o 24_06_UP_H4K20me2_matrix.gz --skipZeros --outFileNameMatrix 24_06_UP_H3K9me2_scaled.tab --outFileSortedRegions 24_06_UP_H4K20me2_sorted.bed
#By default, computeMatrix will handle the flipping for the negative strand, ensuring that TSS is aligned for genes on both strands in your plotProfile or plotHeatmap visualizations.
plotProfile -m 24_06_UP_H4K20me2_matrix.gz -out 24_06_UP_H4K20me2.pdf --plotTitle "H4K20me2 distributions" --plotType=se --yMin=0.8 --yMax=1.8 --perGroup --plotHeight 10 --plotWidth 12 --legendLocation best --plotFileFormat pdf


###### DOWN DEG
## H3K4me2
computeMatrix reference-point -S H3K4me2_WT_IgG_ratio.bw H3K4me2_KO_IgG_ratio.bw -R 2024_06_19_DOWN_DEG.bed -b 30000 -a 30000 --referencePoint TSS -o 24_06_DOWN_H3K4me2_matrix.gz --skipZeros --outFileNameMatrix 24_06_DOWN_H3K4me2_scaled.tab --outFileSortedRegions 24_06_DOWN_H3K4me2_sorted.bed
#By default, computeMatrix will handle the flipping for the negative strand, ensuring that TSS is aligned for genes on both strands in your plotProfile or plotHeatmap visualizations.
plotProfile -m 24_06_DOWN_H3K4me2_matrix.gz -out 24_06_DOWN_H3K4me2.pdf --plotTitle "H3K4me2 distributions" --plotType=se --yMin=0.8 --yMax=1.9 --perGroup --plotHeight 10 --plotWidth 12 --legendLocation best --plotFileFormat pdf
## H3K9me2
computeMatrix reference-point -S H3K9me2_WT_IgG_ratio.bw H3K9me2_KO_IgG_ratio.bw -R 2024_06_19_DOWN_DEG.bed -b 30000 -a 30000 --referencePoint TES -o 24_06_DOWN_H3K9me2_matrix.gz --skipZeros --outFileNameMatrix 24_06_DOWN_H3K9me2_scaled.tab --outFileSortedRegions 24_06_DOWN_H3K9me2_sorted.bed
#By default, computeMatrix will handle the flipping for the negative strand, ensuring that TSS is aligned for genes on both strands in your plotProfile or plotHeatmap visualizations.
plotProfile -m 24_06_DOWN_H3K9me2_matrix.gz -out 24_06_DOWN_H3K9me2.pdf --plotTitle "H3K9me2 distributions" --plotType=se --yMin=0.8 --yMax=1.8 --perGroup --plotHeight 10 --plotWidth 12 --legendLocation best --plotFileFormat pdf
## H4K20me2
computeMatrix reference-point -S H4K20me2_WT_IgG_ratio.bw H4K20me2_KO_IgG_ratio.bw -R 2024_06_19_DOWN_DEG.bed -b 30000 -a 30000 --referencePoint TES -o 24_06_DOWN_H4K20me2_matrix.gz --skipZeros --outFileNameMatrix 24_06_DOWN_H3K9me2_scaled.tab --outFileSortedRegions 24_06_DOWN_H4K20me2_sorted.bed
#By default, computeMatrix will handle the flipping for the negative strand, ensuring that TSS is aligned for genes on both strands in your plotProfile or plotHeatmap visualizations.
plotProfile -m 24_06_DOWN_H4K20me2_matrix.gz -out 24_06_DOWN_H4K20me2.pdf --plotTitle "H4K20me2 distributions" --plotType=se --yMin=0.8 --yMax=1.8 --perGroup --plotHeight 10 --plotWidth 12 --legendLocation best --plotFileFormat pdf


#### Constant
## H3K4me2
computeMatrix reference-point -S H3K4me2_WT_IgG_ratio.bw H3K4me2_KO_IgG_ratio.bw -R constant_non_DEG.bed -b 30000 -a 30000 --referencePoint TSS -o 24_09_Cons_H3K4me2_matrix.gz --skipZeros --outFileNameMatrix 24_09_Cons_H3K4me2_scaled.tab --outFileSortedRegions 24_09_Cons_H3K4me2_sorted.bed
#By default, computeMatrix will handle the flipping for the negative strand, ensuring that TSS is aligned for genes on both strands in your plotProfile or plotHeatmap visualizations.
plotProfile -m 24_09_Cons_H3K4me2_matrix.gz -out 24_09_Cons_H3K4me2.pdf --plotTitle "H3K4me2 distributions" --plotType=se --yMin=0.8 --yMax=1.9 --perGroup --plotHeight 10 --plotWidth 12 --legendLocation best --plotFileFormat pdf
## H3K9me2
computeMatrix reference-point -S H3K9me2_WT_IgG_ratio.bw H3K9me2_KO_IgG_ratio.bw -R constant_non_DEG.bed -b 30000 -a 30000 --referencePoint TES -o 24_09_Cons_H3K9me2_matrix.gz --skipZeros --outFileNameMatrix 24_09_Cons_H3K9me2_scaled.tab --outFileSortedRegions 24_09_Cons_H3K9me2_sorted.bed
#By default, computeMatrix will handle the flipping for the negative strand, ensuring that TSS is aligned for genes on both strands in your plotProfile or plotHeatmap visualizations.
plotProfile -m 24_09_Cons_H3K9me2_matrix.gz -out 24_09_Cons_H3K9me2.pdf --plotTitle "H3K9me2 distributions" --plotType=se --yMin=0.8 --yMax=1.8 --perGroup --plotHeight 10 --plotWidth 12 --legendLocation best --plotFileFormat pdf
## H4K20me2
computeMatrix reference-point -S H4K20me2_WT_IgG_ratio.bw H4K20me2_KO_IgG_ratio.bw -R constant_non_DEG.bed -b 30000 -a 30000 --referencePoint TES -o 24_09_Cons_H4K20me2_matrix.gz --skipZeros --outFileNameMatrix 24_09_Cons_H3K9me2_scaled.tab --outFileSortedRegions 24_09_Cons_H4K20me2_sorted.bed
#By default, computeMatrix will handle the flipping for the negative strand, ensuring that TSS is aligned for genes on both strands in your plotProfile or plotHeatmap visualizations.
plotProfile -m 24_09_Cons_H4K20me2_matrix.gz -out 24_09_Cons_H4K20me2.pdf --plotTitle "H4K20me2 distributions" --plotType=se --yMin=0.8 --yMax=1.8 --perGroup --plotHeight 10 --plotWidth 12 --legendLocation best --plotFileFormat pdf
