#!/bin/bash

#SBATCH --job-name=go_1
#SBATCH --nodes=1
#SBATCH --time=12:00:00
#SBATCH --array=1-916
#SBATCH --mem=10000M
#SBATCH --output=slurm-%A_%a.out
#SBATCH --partition=serial

echo "Running on ${HOSTNAME}"
module add languages/r/3.6.2
Rscript all_network_PheWAS_analysis_iterative_GO.R /mnt/storage/scratch/hw15842/repo/Chapter5_BMI_triplicate_analyses/Data/ /mnt/storage/scratch/hw15842/repo/Chapter5_BMI_triplicate_analyses/Data/triplicate_iterative_results/ --row_num=${SLURM_ARRAY_TASK_ID} 4000