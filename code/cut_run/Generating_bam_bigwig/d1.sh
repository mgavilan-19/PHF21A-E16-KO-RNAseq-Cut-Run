#!/usr/bin/bash

#SBATCH --job-name=cut_run_bt2_bw          ## Name of the job for the scheduler
#SBATCH --account=xxxxxx         ## Your PI's uniqname plus 0,99, or other number
#SBATCH --partition=largemem         ## Choose: standard, largemem, gpu, spgpu, debug            
##SBATCH --gpus=1  			   ## if partition=gpu, number of GPUS needed
#SBATCH --nodes=1                      ## number of nodes you are requesting
#SBATCH --ntasks=1                     ## how many task spaces do you want to reserve
#SBATCH --cpus-per-task=16            ## how many cores do you want to use per task
#SBATCH --time=24:00:00                 ## Maximum length of time you are reserving the 
#SBATCH --mem=120G                      ## Memory limit  
#SBATCH --mail-user=mgavilan@umich.edu  ## send email notifications to umich email listed
#SBATCH --mail-type=ALL                ## when to send email (standard values are:
                                       ## NONE,BEGIN,END,FAIL,REQUEUE,ALL.
   ## (See documentation for others)
#SBATCH --output=./%x-%j               ## send output and error info to the file listed
                                       ##(optional: different name format than default) 

# activate conda enviroment to use snakemakes
source ~/miniconda3/etc/profile.d/conda.sh # change if needed
conda activate /home/mgavilan/miniconda3/envs/snakemake # change

# make nessesary directories
mkdir -p bam_shigeki
mkdir -p bigwig_shigeki

# run snakemake to map the bowtie2 output to bigwig
snakemake --cores 16 -s d1_mapping_to_bw.smk --rerun-incomplete 
