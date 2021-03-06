#!/bin/bash

#SBATCH --job-name=go_1
#SBATCH --nodes=1
#SBATCH --time=12:00:00
#SBATCH --array=1-1000%1
#SBATCH --mem=10000M
#SBATCH --output=slurm-%A_%a.out
#SBATCH --partition=serial

echo "Running on ${HOSTNAME}"
module add languages/r/3.6.2
Rscript triplicate_MR_iterative_script.R /mnt/storage/scratch/hw15842/repo/Chapter5_BMI_triplicate_analyses/Data/ /mnt/storage/scratch/hw15842/repo/Chapter5_BMI_triplicate_analyses/Data/triplicate_iterative_results/  --row_num=${SLURM_ARRAY_TASK_ID} 1000