# !/bin/bash
#SBATCH --job-name=seqtk_samples
#SBATCH --mail-user= xxxx
#SBATCH --mail-type=BEGIN,END
#SBATCH --account=xxxxxx
#SBATCH --partition=largemem
#SBATCH --time=7:00:00
#SBATCH --mem=163840M

source ~/miniconda3/etc/profile.d/conda.sh 
conda activate seqtk-env

# Random sampling

# IgGko2_rs
seqtk sample 9932-BZ-4_S3_R1_001.fastq.gz 30000000 > /nfs/turbo/umms-siwase/CECILIA/new_Cut_Run_Bo/raw_fastq/9932-BZ-4_S3_R1_001_t_30.fastq
seqtk sample 9932-BZ-4_S3_R2_001.fastq.gz 30000000 > /nfs/turbo/umms-siwase/CECILIA/new_Cut_Run_Bo/raw_fastq/9932-BZ-4_S3_R2_001_t_30.fastq

# H3K9me2wt1_rs
seqtk sample 9932-BZ-9_S4_R1_001.fastq.gz 30000000 > /nfs/turbo/umms-siwase/CECILIA/new_Cut_Run_Bo/raw_fastq/9932-BZ-9_S4_R1_001_t_30.fastq
seqtk sample 9932-BZ-9_S4_R2_001.fastq.gz 30000000 > /nfs/turbo/umms-siwase/CECILIA/new_Cut_Run_Bo/raw_fastq/9932-BZ-9_S4_R2_001_t_30.fastq

# H4K20me2ko1
seqtk sample 9932-BZ-15_S5_R1_001.fastq.gz 30000000 > /nfs/turbo/umms-siwase/CECILIA/new_Cut_Run_Bo/raw_fastq/9932-BZ-15_S5_R1_001_t_30.fastq
seqtk sample 9932-BZ-15_S5_R2_001.fastq.gz 30000000 > /nfs/turbo/umms-siwase/CECILIA/new_Cut_Run_Bo/raw_fastq/9932-BZ-15_S5_R2_001_t_30.fastq
