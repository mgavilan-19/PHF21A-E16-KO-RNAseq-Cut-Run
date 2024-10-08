#!/bin/bash
#SBATCH --job-name=average_bw
#SBATCH --mail-user=XXXXXXXX
#SBATCH --mail-type=BEGIN,END
#SBATCH --account=XXXXXXXX
#SBATCH --partition=largemem
#SBATCH --time=3:00:00
#SBATCH --mem=163840M

module load Bioinformatics
module load python 
source ~/miniconda3/etc/profile.d/conda.sh  # change if needed
conda activate deeptools

#IgG - WT
bigwigCompare -b1 IgWT1.bw  -b2 IgWT2.bw  -o IgG_WT_reps_merged.bw --operation mean 
#IgG - KO
bigwigCompare -b1 IgKO1.bw  -b2 IgKO2.bw  -o IgG_KO_reps_merged.bw --operation mean
# H3K4me2 - WT
bigwigCompare -b1 h3k4me2_wt1.bw -b2 h3k4me2_wt2.bw -o H3K4me2_WT_reps_merged.bw --operation mean 
# H3K4me2 - KO
bigwigCompare -b1 h3k4me2_ko1.bw -b2 h3k4me2_ko2.bw -o H3K4me2_KO_reps_merged.bw --operation mean 
# H3K9me2 - WT
bigwigCompare -b1 h3k9me2_wt1.bw -b2 h3k9me2_wt2.bw -o H3K9me2_WT_reps_merged.bw --operation mean 
# H3K9me2 - KO
bigwigCompare -b1 h3k9me2_ko1.bw -b2 h3k9me2_ko2.bw -o H3K9me2_KO_reps_merged.bw --operation mean 
# H4K20me2 - WT
bigwigCompare -b1 h4k20me2_wt1.bw -b2 h4k20me2_wt2.bw -o H4K20me2_WT_reps_merged.bw --operation mean 
# H4K20me2 - KO
bigwigCompare -b1 h4k20me2_ko1.bw -b2 h4k20me2_ko2.bw -o H4K20me2_KO_reps_merged.bw --operation mean 

