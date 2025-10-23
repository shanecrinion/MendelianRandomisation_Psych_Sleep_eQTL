# #!/usr/bin/env Rscript 

# Shane Crinion / shanecrinion@gmail.com / s.crinion1@nuigalway.ie
# 11-2-2022
# This script processes ugly non-standard GWAS summary stats and formats for TwoSampleMR

# PROCESSES & OUTPUTS
# - clean and tidy instrument data to use as exposure
# ├── import list of top loci 
# ├── format data
# output: processed list of instruments 
# - clean and tidy instrument data to use as outcome
# ├── import summary stats
# ├── format data
# output: processed summary stats
# - create metadata file
# ├── collect data name, description, source & significant loci counts
# output: metadata file for datasets used


# import libraries
suppressPackageStartupMessages({
  library(data.table) # import data
  library(TwoSampleMR)})

# set up metadata file
header <- paste("data", "consortia", "pop", "cases","controls","loci.5e-8", sep = ',')
write(header,file="gwas_processed/gwas_exposure_metadata.csv", append=TRUE)

write_metadata <- function(){
  line <- paste(name, corsortia,  pop, cases, controls, loci, source, sep=",")
  write(line,file="gwas_processed/gwas_exposure_metadata.csv", append=TRUE)}


## 1. SCZ 
# info for metadata file
name <- "Schizophrenia"
consortia <- "PGC"
cases <- 53386
controls <- 77258
source <- "Wave 3 EU only subset requested directly from the PGC"
description <- "SNPs reported as GWS in a meta-analysis of EU+EA were used, however statistics were extracted from a EU only subset."

# for sz , GWS loci are extracted from one table and then SE and Beta are extracted from EU-only summary stats
snp_list <- read.csv('PGC3-supp/Supplementary Table 3 - Combined discovery-replication loci.csv')$top.index
exposure_data <- read_exposure_data("PGC3-supp/EU/daner_PGC_SCZ_w3_90_0518d_eur.gz", 
                     sep= '\t', snp_col = "SNP",beta_col = "OR", se_col = "SE", 
                     effect_allele_col = "A1", other_allele_col = "A2", 
                     eaf_col = "FRQ_A_53386",   pval_col = "P")
exposure_data$beta.exposure <- log(exposure_data$beta.exposure) # use log(OR)for binary traits


