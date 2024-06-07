#!/usr/bin/env Rscript
library(data.table)
library(tidyverse)
library(vautils)

setwd("/mnt/sda/gagelo01/Projects/F2_F11/")
ID<-"dis-1-1"
vcffile<-paste0("/mnt/sda/gagelo01/Vcffile/Server_vcf/", ID, "/", ID, ".vcf.gz")
gwasvcf::set_bcftools()
gwasvcf::set_plink()
inst <- GagnonMR::get_inst(vcffile = vcffile, pval = 5e-8, clump = TRUE, r2 = 0.01)
input_data <- data.table(inst[,.(SNP, chr.exposure, pos.exposure)])
setnames(input_data, c("SNP", "chr.exposure", "pos.exposure"), c("rsid", "chromosome", "position"))

genes <- find_nearest_gene(input_data,flanking=100,build="hg19",collapse=FALSE)
setDT(genes)
genes[,distance := ifelse(distance=="intergenic","0",distance)%>%as.numeric(.)]
genes[GENE %in% c("F2", "F11"), ]

genes <- genes[, GENE[which.min(abs(distance))], by = "rsid"]
setnames(genes, "V1", "nearest_gene")      

source("/mnt/sda/gagelo01/Projects/Brain_pQTL/Analysis/CMplot_eloi.R")
# tomanhatthan<- res_map[ , .SD[which.max(posprob_coloc_PPH4)] , by = "hgnc_symbol", .SDcols = c("posprob_colocH4.SNP")]
# data <- VariantAnnotation::readVcf("/mnt/sda/gagelo01/Vcffile/Server_vcf/dis-1-1/dis-1-1.vcf.gz") %>%
#   gwasglue::gwasvcf_to_TwoSampleMR(., type = "outcome") %>% as.data.table(.)
data <- gwasvcf::query_gwas(vcf = vcffile, pval = 5e-5) %>%
  gwasglue::gwasvcf_to_TwoSampleMR(., type = "outcome") %>% as.data.table(.)
data <- data[!(is.na(chr.outcome) | is.na(pos.outcome) | is.na(pval.outcome)),]
data <- data[,.(SNP, chr.outcome, pos.outcome, pval.outcome)]

# setwd(paste0(wd, "/Results"))
# 
# 
# CMplot_eloi(as.data.frame(data), type="p", plot.type="m", LOG10=TRUE, threshold=5e-08, amplify=FALSE, memo="", dpi=300, verbose=TRUE, width=14,
#             height=6, col=c("grey70", "grey90"), threshold.lwd=2, cex=0.4, highlight.cex =  0.6, highlight=genes$rsid,
#             highlight.text=genes$nearest_gene, highlight.col = "#9E131E",
#             highlight.text.xadj =rep(-1, genes[,.N]), highlight.text.yadj =rep(1, genes[,.N]),
#             file.output = TRUE, arrow_col_eloi = "black", highlight.text.cex=0.5, padding = 2, file=c("tiff"))

setwd(paste0(wd))