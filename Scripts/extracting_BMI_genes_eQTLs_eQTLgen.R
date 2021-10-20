################################################
### extracting_BMI_genes_eQTLs_eQTLgen.R ###
################################################

## module add languages/r/4.0.3 
# has to be in R version 4

args  <-  commandArgs(trailingOnly=TRUE)
data_location <- toString(args[1])

# data_location <- "/mnt/storage/scratch/hw15842/repo/Chapter5_BMI_triplicate_analyses/Data/"

library("data.table")
library("dplyr")
library("plyr")


setwd(paste0(data_location))



sig_trans_eQTLs_filename <- gzfile('2018-09-04-trans-eQTLsFDR0.05-CohortInfoRemoved-BonferroniAdded.txt.gz','rt')  
sig_trans_eQTLs <- read.table(sig_trans_eQTLs_filename,header=T)


sig_cis_eQTLs_filename <- gzfile('2019-12-11-cis-eQTLsFDR0.05-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt.gz','rt')  
sig_cis_eQTLs <- read.table(sig_cis_eQTLs_filename,header=T)


all_sig_eQTLs <- rbind(sig_cis_eQTLs, sig_trans_eQTLs)


load("BMI_gene_list.rdata")

BMI_gene_eQTLs <- subset(all_sig_eQTLs, all_sig_eQTLs$GeneSymbol %in% BMI_gene_list)


save(BMI_gene_eQTLs, file="BMI_gene_eQTLs.rdata")