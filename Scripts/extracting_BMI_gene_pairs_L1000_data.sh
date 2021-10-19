#!/bin/bash

#SBATCH --job-name=l1000
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=20
#SBATCH --time=48:00:00
#SBATCH --mem=100000M
#SBATCH --partition=cpu


echo "Running on ${HOSTNAME}"
module add languages/r/4.0.3 
Rscript extracting_BMI_gene_pairs_L1000_data.R /mnt/storage/scratch/hw15842/repo/Chapter5_BMI_triplicate_analyses/Data/