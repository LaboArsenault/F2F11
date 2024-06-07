#!/usr/bin/env Rscript
setwd("/mnt/sda/gagelo01/Projects/F2_F11")
source("Analysis/tosource.R")

dt_gene_region <- fread("/mnt/sda/gagelo01/Projects/small_MR_exploration/Test_potential_project/Tribal/Data/Modified/dt_gene_region.txt")
coagulation_gene <- fread("Data/Raw/bloodcoagulationpathway_genes.txt") #genes from https://maayanlab.cloud/Harmonizome/gene_set/Blood+coagulation/PANTHER+Pathways
gene_toinclude <-  dt_gene_region[hgnc %in% coagulation_gene$Symbol, hgnc %>% unique]

dt_gene_region <- dt_gene_region[hgnc%in%gene_toinclude,]
dt_gene_region[,gene_region:=NULL]
dt_gene_region[, gene_region := paste0(chr, ":", (start-1e6)%>%ifelse(.<1, 1, .), "-", (end+1e6))]

fwrite("Data/Modified/dt_gene_region.txt")
message("This script finished without errors")

# setwd("/mnt/sda/gagelo01/Projects/F2_F11")
# source("Analysis/tosource.R")
# dt_gene_region <- fread("/mnt/sda/gagelo01/Projects/small_MR_exploration/Test_potential_project/PWMR_anything/Data/Modified/dt_gene_region.txt")
# dt_gene_region[study=="Atherosclerosis Risk in Communities", study := "ARIC"]
# 
# k<- df_index[grepl("eqtl-4-", id), .(id, clean_variable_name, consortium, trait)]
# setnames(k, c("clean_variable_name", "consortium"), c("hgnc","study"))
# k<- merge(k, distinct(dt_gene_region[,.(hgnc, UniProt, gene_region)]), by = "hgnc")
# dt_gene_region <- rbind(dt_gene_region, k)
# dt_gene_region[grepl("eqtl-74", id), vcffile := paste0(wd, "/Data/Modified/Gtex_vcf/", id, "/", id, ".vcf.gz")]
# dt_gene_region[!grepl("eqtl-74", id), vcffile := paste0("/mnt/sda/gagelo01/Vcffile/Server_vcf/", id, "/", id, ".vcf.gz")]
# 
# fwrite(dt_gene_region, "Data/Modified/dt_gene_region.txt")
# message("This script finished without errors")