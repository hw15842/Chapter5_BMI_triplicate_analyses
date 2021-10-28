######################################################
#### triplicate_MR_iterative_script.R ####
######################################################

args  <-  commandArgs(trailingOnly=TRUE)
data_location <- toString(args[1])
results_location <- toString(args[2])
row_num <- toString(args[3])
batch <- toString(args[4])


# data_location <- "/mnt/storage/scratch/hw15842/repo/Chapter5_BMI_triplicate_analyses/Data/"
# results_location<- "/mnt/storage/scratch/hw15842/repo/Chapter5_BMI_triplicate_analyses/Data/triplicate_iterative_results/"


library(TwoSampleMR)
library(gwasglue)
library(gwasvcf)
library(ieugwasr)
library(dplyr)
library(plyr)
library(data.table)

load(paste0(data_location,"sig_l1000_BMI_gene_pairs_to_use.rdata"))

row <- unlist(regmatches(row_num, gregexpr('\\(?[0-9,.]+', network_num))) ### This stops it pasting the full "--section_of_subset=1" and just keeps the "1" 

row <- as.numeric(noquote(row)) + as.numeric(batch)

row 




MR_func <- function(gene_pair_line){
  
  df <- sig_l1000_BMI_gene_pairs[gene_pair_line,]
  print(df)
  gene1_ENSG_ID <- subset(BMI_genes_ENSG_ID, BMI_genes_ENSG_ID$GeneSymbol == df$gene1)[1,]
  print(gene1_ENSG_ID)
  gene2_ENSG_ID <- subset(BMI_genes_ENSG_ID, BMI_genes_ENSG_ID$GeneSymbol == df$gene2)[1,]
  
  # Need to extract them from ieugwasr so can harmonise
  gene1_exp <- extract_instruments(paste0("eqtl-a-",gene1_ENSG_ID$Gene))
  gene2_exp <- extract_instruments(paste0("eqtl-a-",gene2_ENSG_ID$Gene))
  
  
  # remove any matching SNPs
  snps1 <- gene1_exp$SNP
  snps2 <- gene2_exp$SNP
  gene1_exp <- subset(gene1_exp, !gene1_exp$SNP %in% snps2)
  gene2_exp <- subset(gene2_exp, !gene2_exp$SNP %in% snps1)

  
  # extract snps from BMI
  gene1_out <- extract_outcome_data(gene1_exp$SNP, "ukb-a-248", proxies = F)
  gene2_out <- extract_outcome_data(gene2_exp$SNP, "ukb-a-248", proxies = F)
  gene1_on_gene_2_out <- extract_outcome_data(gene1_exp$SNP, paste0("eqtl-a-",gene2_ENSG_ID$Gene), proxies=F)

  
  # harmonise the data if both have the data frame 
  if(is.null(gene1_out) == "FALSE" & 
     is.null(gene2_out) == "FALSE" & 
     is.null(gene1_on_gene_2_out) == "FALSE"){
  gene1_dat <- harmonise_data(gene1_exp, gene1_out)
  gene2_dat <- harmonise_data(gene2_exp, gene2_out)
  g1_on_g2_dat <- harmonise_data(gene1_exp, gene1_on_gene_2_out)
  
  # run MR
  gene1_MR <- mr(gene1_dat)
  gene2_MR <- mr(gene2_dat)
  g1_on_g2_MR <- mr(g1_on_g2_dat)
  
  colnames(gene1_MR) <- paste("Gene1", colnames(gene1_MR), sep = "_")
  colnames(gene2_MR) <- paste("Gene2", colnames(gene2_MR), sep = "_")
  colnames(g1_on_g2_MR) <- paste("G1_on_G2", colnames(g1_on_g2_MR), sep = "_")
  
  res <- cbind(gene1_MR, gene2_MR, g1_on_g2_MR)
  res$gene1 <- df$gene_1
  res$gene1_chr <- gene1_ENSG_ID$GeneChr
  res$gene1_pos <- gene1_ENSG_ID$GenePos
  
  res$gene2 <- df$gene_2
  res$gene2_chr <- gene2_ENSG_ID$GeneChr
  res$gene2_pos <- gene2_ENSG_ID$GenePos
  }
  else{ 
    res <- df
    res$gene1_chr <- gene1_ENSG_ID$GeneChr
    res$gene1_pos <- gene1_ENSG_ID$GenePos
    res$gene2_chr <- gene2_ENSG_ID$GeneChr
    res$gene2_pos <- gene2_ENSG_ID$GenePos
    res <- res %>% add_row() 
    res <- res %>% add_row() 
    res <- res %>% add_row() 
    res <- res %>% add_row() 
    res <- res %>% add_row() 
    res <- res %>% add_row() 
    # adding in 6 rows of NAs so that can search for the dataframes with 7 rows as the ones that did not have the data to run
    # dataframes with 1 row are where there is only the wald method applied 
    # data frames with 5 rows are ones where there is an IVW or the other methods as well
    # worried that there may be an occasion where if all of them IVW might have 2 rows so making 7 rows just to be certain 
  }

  return(res)

}



MR_res <- MR_func(row)

save(MR_res, file=paste0(results_location, "MR_results_", row))



