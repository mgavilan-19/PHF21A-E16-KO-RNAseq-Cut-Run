# !/bin/bash
#SBATCH --job-name=bigwig_ratio
#SBATCH --mail-user=xxxxx
#SBATCH --mail-type=BEGIN,END
#SBATCH --account=xxxxxx
#SBATCH --partition=largemem
#SBATCH --time=5:00:00
#SBATCH --mem=163840M

module load python 
source ~/miniconda3/etc/profile.d/conda.sh  # change if needed
conda activate deeptools

#This tool compares two bigWig files: treated/control samples

# H3K4me2 and IgG ---> WT
bigwigCompare -b1 H3K4me2_WT_reps_merged.bw -b2  IgG_WT_reps_merged.bw --operation=ratio -o H3K4me2_WT_IgG_ratio.bw --pseudocount 1

# H3K4me2 and IgG ---> KO
bigwigCompare -b1 H3K4me2_KO_reps_merged.bw -b2 IgG_KO_reps_merged.bw --operation=ratio -o H3K4me2_KO_IgG_ratio.bw --pseudocount 1

# H3K9me2 and IgG ---> WT
bigwigCompare -b1 H3K9me2_WT_reps_merged.bw -b2 IgG_WT_reps_merged.bw --operation=ratio -o H3K9me2_WT_IgG_ratio.bw --pseudocount 1

# H3K9me2 and IgG ---> KO
bigwigCompare -b1 H3K9me2_KO_reps_merged.bw -b2 IgG_KO_reps_merged.bw --operation=ratio -o H3K9me2_KO_IgG_ratio.bw --pseudocount 1

# H4K20me2 and IgG ---> WT
bigwigCompare -b1 H4K20me2_WT_reps_merged.bw -b2 IgG_WT_reps_merged.bw --operation=ratio -o H4K20me2_WT_IgG_ratio.bw --pseudocount 1

# H4K20me2 and IgG ---> KO
bigwigCompare -b1 H4K20me2_KO_reps_merged.bw -b2 IgG_KO_reps_merged.bw --operation=ratio -o H4K20me2_KO_IgG_ratio.bw --pseudocount 1


