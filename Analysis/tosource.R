#!/usr/bin/env Rscript
library(data.table)
library(tidyverse)
library(GagnonMR)
library(tictoc)
library(furrr)

wd<-"/mnt/sda/gagelo01/Projects/F2_F11"
setwd(wd)
gwasvcf::set_bcftools()
gwasvcf::set_plink()
ldref<-"/home/couchr02/Mendel_Commun/Christian/LDlocal/EUR_rs"
x <- paste0("awk -F ' ' '{print $2}' ","/home/couchr02/Mendel_Commun/Christian/LDlocal/EUR_rs", ".bim")
dt_ref_snp <- fread(cmd = x, header = FALSE) #Those are SNP in 1000G with MAF>0.01
df_index <- data.table::fread("/mnt/sda/gagelo01/Vcffile/server_gwas_id.txt")
ao <- fread("/mnt/sda/gagelo01/Vcffile/available_outcomes_2021-10-13.txt")
ao_small <- ao[id %in% list.files("/mnt/sda/gagelo01/Vcffile/MRBase_vcf/"), ]
dt_gene_region <- fread("/mnt/sda/gagelo01/Projects/small_MR_exploration/Test_potential_project/PWMR_anything/Data/Modified/dt_gene_region.txt")

options(future.globals.maxSize= 5e9)
plan(multicore, workers = 40, gc = TRUE)

###chose instrument
study_to_selectinst <- list("deCODE", "FENLAND", "ARIC", "IUCPQ Biobank")
holistic_selection<-FALSE #if TRUE will use information on all study. If FALSE will not.
study_noinst_butexposure <- NULL
should_skip_homogenous = TRUE
typeof_sel = c("lead_snp", "multicis_independent_clumping")
should_select_pan <- TRUE
######change parameters#########
######
parameters <- GagnonMR::default_param()
parameters$path <- c(GagnonMR::default_param()$path, paste0(wd, "/Data/Modified/Gtex_vcf/"))
parameters$uni_cis_minpval <- 1
parameters$snp_bim <- dt_ref_snp$V1 #I only include SNP in 1000G with maf > 0.01
parameters$multicis_clumping["clump_r2"]<-0.4
parameters$proxy_rsq <- 0.6
############Choose the gene you wish to include and the outcome you wish to include  ###########
coagulation_gene <- fread("Data/Raw/bloodcoagulationpathway_genes.txt") #genes from https://maayanlab.cloud/Harmonizome/gene_set/Blood+coagulation/PANTHER+Pathways
gene_toinclude <-  dt_gene_region[hgnc %in% coagulation_gene$Symbol, hgnc %>% unique]
vec_tissue_gtex <- "Liver" 
############Choose the outcome you wish to include  ###########
ID_server_out <- c(
  "dis-15-157", "dis-15-1516", "dis-15-1739",#various bleeding outcome with ncase < 1000 in Finngen
  "dis-1-1", "trait-7-1",  #venout_thromboembolism, parental lifespant (years),
  df_index[grepl("dis-14-", id) & population == "European",id]) #Ischemicstroke, Cardioembolic stroke
out_server <- paste0("/mnt/sda/gagelo01/Vcffile/Server_vcf/", ID_server_out, "/", ID_server_out, ".vcf.gz")
out_mrbase<-NULL

split_outcome<-FALSE
all_mr_methods_short = FALSE
all_mr_methods_skip_presso = FALSE
tsmr_robust = c( "mr_weighted_median", "mr_egger_regression")