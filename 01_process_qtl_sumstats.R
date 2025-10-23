# #!/usr/bin/env Rscript 

# Shane Crinion / shanecrinion@gmail.com / s.crinion1@nuigalway.ie
# 10-2-2022
# This script cleans non-standardised QTL summary stats and formats for TwoSampleMR

# PROCESSES & OUTPUTS
# - clean and tidy QTL sumstats
# ├── import summary stats, lookup and header files
# ├── extract significant SNPs
# ├── extract annotations from lookup file
# ├── format data
# ├── clump data
# output: processed sumstats for QTL data
# - create metadata file
# ├── collect data name, description, source & SNP and gene counts
# output: metadata file for datasets used

# DATASETS
# -BRAIN
# - PSYCHencode ✓
# ├── Cortex ✓
# - GTEx Release 7
# ├── Cortex ✓
# ├── Cerebellum
# ├── Hypothalamus ✓ 
# ├── Whole Blood  ✓


# Libraries
suppressPackageStartupMessages({
   library(data.table) 
   library(TwoSampleMR)})

# Create object of meta-data 
header <- paste("data", "region", "individuals","individuals.descr",
                "count.variants","count.5e-8","count.duplicated",
                "snp.gene.pairs","snps.unique",
                "genes.unique", "source", "data.description", sep = ',')
write(header,file="metadata.csv", append=TRUE)

# Write function to write meta-data
write_metadata <- function(){
   line <- paste(name, region, individuals,  individuals_descr, 
                 count.variants, count.sig, count.duplicated, 
                 count.snps, count.snps.unique, 
                 count.genes.unique, source, data_descr, sep=",")
   write(line,file="metadata.csv", append=TRUE)}

## 1. PSYCHencode 
# Info for metadata file
name <- "PSYCHencode"
region <- "prefrontal cortex"
individuals <- 1387
individuals_descr <- "Filtered adult samples with matching gene expression and genotypes (679 healthy controls + 497 schizophrenia + 172 bipolar disorder + 31 autism spectrum disorder and 8 affective disorder patients)" 
data_descr <-  "Full set of cis-eQTLs with no p-value or FDR filtering"
source <- "http://resource.psychencode.org/Datasets/Derived/QTLs/Full_hg19_cis-QTL.txt.gz"

# Import data
QTL <- fread('psychencode/psychencode_Full_hg19_cis-QTL.txt.gz',sep = ' ')
colnames(QTL) <- names(fread("psychencode/DER-08b_hg19_QTL.bonferroni.txt", sep="\t"))[1:14]
QTL.lookup <- fread('psychencode/SNP_Information_Table_with_Alleles.txt', sep="\t")
QTL.sig <- subset(QTL, nominal_pval < 5e-8)

# get RSID and alleles for remaining SNPs
QTL.sig$rsid <- QTL.lookup$Rsid[match(QTL.sig$SNP_id, QTL.lookup$PEC_id)] 
QTL.sig$ref <- QTL.lookup$REF[match(QTL.sig$SNP_id, QTL.lookup$PEC_id)]
QTL.sig$alt <- QTL.lookup$ALT[match(QTL.sig$SNP_id, QTL.lookup$PEC_id)]

# calculate SE
tail <- 2
se <- abs(QTL.sig$regression_slope/ qnorm(QTL.sig$nominal_pval/tail))
QTL.sig$se <- se

# save counts for metadata file
count.variants <- dim(QTL)[1] 
count.sig <- dim(QTL.sig)[1]
count.duplicated <- sum(duplicated(QTL.sig$rsid))

# create a unique locator for each gene-SNP pair (many SNPs associated with > 1 gene)
QTL.sig$locator <- paste(QTL.sig$gene_id, QTL.sig$rsid, QTL.sig$SNP_id, sep="-")


# format for TwoSampleMR (locator col used for SNPs to preserve all unique gene-SNP pairs)
exposure_data <- format_data(dat = QTL.sig, type = "exposure",snp_col ="locator",
                             beta_col = "regression_slope",se_col =  "se",
                             effect_allele_col = "alt",other_allele_col = "ref",pval_col = "nominal_pval")

# for clumping, set "pos" as locator and SNP as "rsid" (hacky but package doesn't preserve custom columns)
exposure_data$pos <- exposure_data$SNP
exposure_data$SNP <- unlist(lapply(strsplit(exposure_data$pos, "-"), '[[' ,2))

# save processed data
write.csv(exposure_data, file = 'processed_PSYCHencode_QTL.csv')
rm(QTL, QTL.lookup, QTL.sig, se)

# clump data 
exposure_data <- clump_data(dat = exposure_data, clump_kb = 10000, clump_r2 = 0.001)
write.csv(exposure_data, file = 'processed_PSYCHencode_QTL.clumped.csv')

# get gene info
exposure_data$entrez <- toupper(unlist(lapply(strsplit(exposure_data$pos, "-"), '[[' ,1))) 
exposure_data$entrez <- unlist(lapply(strsplit(exposure_data$entrez, "[.]"), '[[', 1))

# snp and gene counts
count.snps <- length(exposure_data$SNP) 
count.snps.unique <- length(unique(exposure_data$SNP))
count.genes.unique <- length(unique(exposure_data$entrez))

# write metadata to file
write_metadata()

# 2. GTEx 
# CORTEX
# info for metadata file
name <- "GTEx Analysis V7"
region <- "Frontal cortex"
individuals <- 0
individuals_descr <- "TBC" 
data_descr <-  "TBC"
source <- "https://gtexportal.org/home/datasets"

# import data
QTL <- fread('gtex-brain/Brain_Cortex.allpairs.txt.gz',sep = '\t')
QTL.lookup <- fread('gtex-brain/GTEx_Analysis_2016-01-15_v7_WholeGenomeSeq_635Ind_PASS_AB02_GQ20_HETX_MISS15_PLINKQC.lookup_table.txt.gz', sep="\t")
QTL.sig <- subset(QTL, pval_nominal < 5e-8)

# get RSID and alleles for remaining SNPs
QTL.sig$rsid <- QTL.lookup$rs_id_dbSNP147_GRCh37p13[match(QTL.sig$variant_id, QTL.lookup$variant_id)] 
QTL.sig$ref <- QTL.lookup$ref[match(QTL.sig$variant_id, QTL.lookup$variant_id)]
QTL.sig$alt <- QTL.lookup$alt[match(QTL.sig$variant_id, QTL.lookup$variant_id)]
QTL.sig$locator <- paste(QTL.sig$gene_id, QTL.sig$rsid, QTL.sig$variant_id, sep="-")

# save counts for metadata file
count.variants <- dim(QTL)[1] 
count.sig <- dim(QTL.sig)[1]
count.duplicated <- sum(duplicated(QTL.sig$rsid))

exposure_data <- format_data(dat = QTL.sig, type = "exposure",
                             snp_col = "locator",beta_col = "slope", 
                             se_col = "slope_se",eaf_col = "maf", 
                             effect_allele_col = "alt",other_allele_col = "ref",
                             pos_col = "locator",pval_col = "pval_nominal")

# for clumping, set "pos" as locator and SNP as "rsid" (hacky but package doesn't preserve custom columns)
exposure_data$pos <- exposure_data$SNP
exposure_data$SNP <- unlist(lapply(strsplit(exposure_data$pos, "-"), '[[' ,2))

# save processed data
write.csv(exposure_data, file = 'processed_GTExCortex_eQTL.csv')
rm(QTL, QTL.lookup, QTL.sig, se)

# clump data 
exposure_data <- clump_data(dat = exposure_data, clump_kb = 10000, clump_r2 = 0.001)
write.csv(exposure_data, file = 'processed_GTExCortex_eQTL.clumped.csv')

# get gene info
exposure_data$entrez <- toupper(unlist(lapply(strsplit(exposure_data$pos, "-"), '[[' ,1))) 
exposure_data$entrez <- unlist(lapply(strsplit(exposure_data$entrez, "[.]"), '[[', 1))

# snp and gene counts
count.snps <- length(exposure_data$SNP) 
count.snps.unique <- length(unique(exposure_data$SNP))
count.genes.unique <- length(unique(exposure_data$entrez))

rm(QTL, QTL.lookup, QTL.sig)

write_metadata()

# info for metadata file
name <- "GTEx Analysis V7"
region <- "Hypothalamus"
individuals <- 0
individuals_descr <- "TBC" 
data_descr <-  "TBC"
source <- "https://gtexportal.org/home/datasets"

# import data
QTL <- fread('gtex-brain/Brain_Hypothalamus.allpairs.txt.gz',sep = '\t')
QTL.lookup <- fread('gtex-brain/GTEx_Analysis_2016-01-15_v7_WholeGenomeSeq_635Ind_PASS_AB02_GQ20_HETX_MISS15_PLINKQC.lookup_table.txt.gz', sep="\t")
QTL.sig <- subset(QTL, pval_nominal < 5e-8)

# get RSID and alleles for remaining SNPs
QTL.sig$rsid <- QTL.lookup$rs_id_dbSNP147_GRCh37p13[match(QTL.sig$variant_id, QTL.lookup$variant_id)] 
QTL.sig$ref <- QTL.lookup$ref[match(QTL.sig$variant_id, QTL.lookup$variant_id)]
QTL.sig$alt <- QTL.lookup$alt[match(QTL.sig$variant_id, QTL.lookup$variant_id)]
QTL.sig$locator <- paste(QTL.sig$gene_id, QTL.sig$rsid, QTL.sig$variant_id, sep="-")

# save counts for metadata file
count.variants <- dim(QTL)[1] 
count.sig <- dim(QTL.sig)[1]
count.duplicated <- sum(duplicated(QTL.sig$rsid))

exposure_data <- format_data(dat = QTL.sig, type = "exposure",
                             snp_col = "locator",beta_col = "slope", 
                             se_col = "slope_se",eaf_col = "maf", 
                             effect_allele_col = "alt",other_allele_col = "ref",
                             pos_col = "locator",pval_col = "pval_nominal")

# for clumping, set "pos" as locator and SNP as "rsid" (hacky but package doesn't preserve custom columns)
exposure_data$pos <- exposure_data$SNP
exposure_data$SNP <- unlist(lapply(strsplit(exposure_data$pos, "-"), '[[' ,2))

# save processed data
write.csv(exposure_data, file = 'processed_GTExHypothalamus_eQTL.csv')
rm(QTL, QTL.lookup, QTL.sig, se)

# clump data 
exposure_data <- clump_data(dat = exposure_data, clump_kb = 10000, clump_r2 = 0.001)
write.csv(exposure_data, file = 'processed_GTExHypthalamus_eQTL.clumped.csv')

# get gene info
exposure_data$entrez <- toupper(unlist(lapply(strsplit(exposure_data$pos, "-"), '[[' ,1))) 
exposure_data$entrez <- unlist(lapply(strsplit(exposure_data$entrez, "[.]"), '[[', 1))

# snp and gene counts
count.snps <- length(exposure_data$SNP) 
count.snps.unique <- length(unique(exposure_data$SNP))
count.genes.unique <- length(unique(exposure_data$entrez))

rm(QTL, QTL.lookup, QTL.sig)

write_metadata()

### updated GTEx data from version 8

# CORTEX
# info for metadata file
name <- "GTEx Analysis V8"
region <- "Cortex"
individuals <- 205
individuals_descr <- "TBC" 
data_descr <-  "Right cerebral frontal pole cortex (sampled at donor collection site and preserved in PAXgene fixative)."
source <- "https://gtexportal.org/home/datasets"

# import data
QTL <- fread('gtexv8/Brain_Cortex.allpairs.txt.gz',sep = '\t')
QTL.sig <- subset(QTL, pval_nominal < 5e-8)
rm(QTL)
QTL.lookup <- fread('gtexv8/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.lookup_table.txt.gz', sep="\t")
head(QTL.lookup)

# get RSID and alleles for remaining SNPs
QTL.sig$rsid <- QTL.lookup$rs_id_dbSNP151_GRCh38p7[match(QTL.sig$variant_id, QTL.lookup$variant_id)] 
QTL.sig$ref <- QTL.lookup$ref[match(QTL.sig$variant_id, QTL.lookup$variant_id)]
QTL.sig$alt <- QTL.lookup$alt[match(QTL.sig$variant_id, QTL.lookup$variant_id)]
QTL.sig$locator <- paste(QTL.sig$gene_id, QTL.sig$rsid, QTL.sig$variant_id, sep="-")
QTL.sig$locator <- gsub("_b38", "", QTL.sig$locator)

# save counts for metadata file
count.variants <- dim(QTL)[1] 
# rm'd QTL so used this - count.variants <- 182641182
count.sig <- dim(QTL.sig)[1]
count.duplicated <- sum(duplicated(QTL.sig$rsid))

exposure_data <- format_data(dat = QTL.sig, type = "exposure",
                             snp_col = "locator",beta_col = "slope", 
                             se_col = "slope_se",eaf_col = "maf", 
                             effect_allele_col = "alt",other_allele_col = "ref",
                             pos_col = "locator",pval_col = "pval_nominal")

# for clumping, set "pos" as locator and SNP as "rsid" (hacky but package doesn't preserve custom columns)
exposure_data$pos <- exposure_data$SNP
exposure_data$SNP <- unlist(lapply(strsplit(exposure_data$pos, "-"), '[[' ,2))

# save processed data
write.csv(exposure_data, file = 'processed_GTExCortexv8_eQTL.csv')
rm(QTL, QTL.lookup, QTL.sig, se)

# clump data 
exposure_data <- clump_data(dat = exposure_data, clump_kb = 10000, clump_r2 = 0.001)
write.csv(exposure_data, file = 'processed_GTExCortexv8_eQTL.clumped.csv')

# get gene info
exposure_data$entrez <- toupper(unlist(lapply(strsplit(exposure_data$pos, "-"), '[[' ,1))) 
exposure_data$entrez <- unlist(lapply(strsplit(exposure_data$entrez, "[.]"), '[[', 1))

# snp and gene counts
count.snps <- length(exposure_data$SNP) 
count.snps.unique <- length(unique(exposure_data$SNP))
count.genes.unique <- length(unique(exposure_data$entrez))

rm(QTL, QTL.lookup, QTL.sig)

write_metadata()

rm(count.duplicated, count.genes.unique, count.sig, 
   count.snps, count.snps.unique, count.variants, 
   data_descr, individuals, individuals_descr, 
   name, region, source)
rm(exposure_data)

# Hippocampus
# info for metadata file
name <- "GTEx Analysis V8"
region <- "Hippocampus"
individuals <- 165
individuals_descr <- "TBC" 
data_descr <-  "Hippocampus (sampled at Miami Brain Bank and preserved as fresh frozen tissue)."
source <- "https://gtexportal.org/home/datasets"

# import data
QTL <- fread('~/Desktop/files/data/gtexv8/Brain_Hippocampus.allpairs.txt.gz',sep = '\t')
# save counts for metadata file
count.variants <- dim(QTL)[1] 
QTL.sig <- subset(QTL, pval_nominal < 5e-8)
rm(QTL)
#write.csv(QTL.sig, file = 'Brain_Hippocampus.allpairs.p5e-8.csv')
#.rs.restartR() # remove large QTL file from memory

#QTL.sig <- read.csv('Brain_Hippocampus.allpairs.p5e-8.csv')
#QTL.lookup <- fread('GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.lookup_table.txt.gz', sep="\t")

# get RSID and alleles for remaining SNPs
QTL.sig$rsid <- QTL.lookup$rs_id_dbSNP151_GRCh38p7[match(QTL.sig$variant_id, QTL.lookup$variant_id)] 
QTL.sig$ref <- QTL.lookup$ref[match(QTL.sig$variant_id, QTL.lookup$variant_id)]
QTL.sig$alt <- QTL.lookup$alt[match(QTL.sig$variant_id, QTL.lookup$variant_id)]
QTL.sig$locator <- paste(QTL.sig$gene_id, QTL.sig$rsid, QTL.sig$variant_id, sep="-")
QTL.sig$locator <- gsub("_b38", "", QTL.sig$locator)

# rm'd QTL so used this - count.variants <- 182641182
count.sig <- dim(QTL.sig)[1]
count.duplicated <- sum(duplicated(QTL.sig$rsid))

exposure_data <- format_data(dat = QTL.sig, type = "exposure",
                             snp_col = "locator",beta_col = "slope", 
                             se_col = "slope_se",eaf_col = "maf", 
                             effect_allele_col = "alt",other_allele_col = "ref",
                             pos_col = "locator",pval_col = "pval_nominal")

# for clumping, set "pos" as locator and SNP as "rsid" (hacky but package doesn't preserve custom columns)
exposure_data$pos <- exposure_data$SNP
exposure_data$SNP <- unlist(lapply(strsplit(exposure_data$pos, "-"), '[[' ,2))

# save processed data
write.csv(exposure_data, file = 'processed_GTExHippocampusv8_eQTL.csv')
rm(QTL, QTL.lookup, QTL.sig, se)

# clump data 
exposure_data <- clump_data(dat = exposure_data, clump_kb = 10000, clump_r2 = 0.001)
write.csv(exposure_data, file = 'processed_GTExHippocampusv8_eQTL.clumped.csv')

# get gene info
exposure_data$entrez <- toupper(unlist(lapply(strsplit(exposure_data$pos, "-"), '[[' ,1))) 
exposure_data$entrez <- unlist(lapply(strsplit(exposure_data$entrez, "[.]"), '[[', 1))

# snp and gene counts
count.snps <- length(exposure_data$SNP) 
count.snps.unique <- length(unique(exposure_data$SNP))
count.genes.unique <- length(unique(exposure_data$entrez))

rm(exposure_data)
rm(QTL.sig)
rm(count.duplicated, count.genes.unique, count.sig, count.snps, count.snps.unique, count.variants, 
   data_descr, individuals, individuals_descr, name,region, source)
rm(exposure_data)

write_metadata()

## hypothalamus
# info for metadata file
name <- "GTEx Analysis V8"
region <- "Hypothalamus"
individuals <- 170
individuals_descr <- "TBC" 
data_descr <-  "Hypothalamus (sampled at Miami Brain Bank and preserved as fresh frozen tissue)."
source <- "https://gtexportal.org/home/datasets"

# import data
QTL <- fread('gtexv8/Brain_Hypothalamus.allpairs.txt.gz',sep = '\t')
# save counts for metadata file
count.variants <- dim(QTL)[1] 
QTL.sig <- subset(QTL, pval_nominal < 5e-8)
rm(QTL)
#.rs.restartR() # remove large QTL file from memory
#write.csv(QTL.sig, file = 'Brain_Hippocampus.allpairs.p5e-8.csv')
# restart R here
#QTL.sig <- read.csv('Brain_Hippocampus.allpairs.p5e-8.csv')

QTL.lookup <- fread('GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.lookup_table.txt.gz', sep="\t")
head(QTL.lookup)

# get RSID and alleles for remaining SNPs
QTL.sig$rsid <- QTL.lookup$rs_id_dbSNP151_GRCh38p7[match(QTL.sig$variant_id, QTL.lookup$variant_id)] 
QTL.sig$ref <- QTL.lookup$ref[match(QTL.sig$variant_id, QTL.lookup$variant_id)]
QTL.sig$alt <- QTL.lookup$alt[match(QTL.sig$variant_id, QTL.lookup$variant_id)]
QTL.sig$locator <- paste(QTL.sig$gene_id, QTL.sig$rsid, QTL.sig$variant_id, sep="-")
QTL.sig$locator <- gsub("_b38", "", QTL.sig$locator)

# rm'd QTL so used this - count.variants <- 182641182
count.sig <- dim(QTL.sig)[1]
count.duplicated <- sum(duplicated(QTL.sig$rsid))

exposure_data <- format_data(dat = QTL.sig, type = "exposure",
                             snp_col = "locator",beta_col = "slope", 
                             se_col = "slope_se",eaf_col = "maf", 
                             effect_allele_col = "alt",other_allele_col = "ref",
                             pos_col = "locator",pval_col = "pval_nominal")

# for clumping, set "pos" as locator and SNP as "rsid" (hacky but package doesn't preserve custom columns)
exposure_data$pos <- exposure_data$SNP
exposure_data$SNP <- unlist(lapply(strsplit(exposure_data$pos, "-"), '[[' ,2))

# save processed data
write.csv(exposure_data, file = 'processed_GTExHypothalamusv8_eQTL.csv')
rm(QTL, QTL.sig, se)

# clump data 
exposure_data <- clump_data(dat = exposure_data, clump_kb = 10000, clump_r2 = 0.001)
write.csv(exposure_data, file = 'processed_GTExHypothalamusv8_eQTL.clumped.csv')

clump_data()

# get gene info
exposure_data$entrez <- toupper(unlist(lapply(strsplit(exposure_data$pos, "-"), '[[' ,1))) 
exposure_data$entrez <- unlist(lapply(strsplit(exposure_data$entrez, "[.]"), '[[', 1))

# snp and gene counts
count.snps <- length(exposure_data$SNP) 
count.snps.unique <- length(unique(exposure_data$SNP))
count.genes.unique <- length(unique(exposure_data$entrez))

rm(exposure_data)
rm(QTL.sig)
rm(count.duplicated, count.genes.unique, count.sig, count.snps, count.snps.unique, count.variants, 
   data_descr, individuals, individuals_descr, name,region, source)
rm(exposure_data)

write_metadata()

# blood

# info for metadata file
name <- "GTEx Analysis V8"
region <- "Whole Blood"
individuals <- 670
individuals_descr <- "TBC" 
data_descr <-  "Femoral vein; subclavian vein and heart are other possible sites."
source <- "https://gtexportal.org/home/datasets"

# import data
library(data.table)
QTL <- fread('~/Desktop/files/data/gtexv8/Whole_Blood.allpairs.txt.gz',sep = '\t')

# save counts for metadata file
count.variants <- dim(QTL)[1] 
QTL.sig <- subset(QTL, pval_nominal < 5e-8)
rm(QTL)
#.rs.restartR() # remove large QTL file from memory
#write.csv(QTL.sig, file = 'Brain_Hippocampus.allpairs.p5e-8.csv')

#QTL.sig <- read.csv('Brain_Hippocampus.allpairs.p5e-8.csv')

QTL.lookup <- fread('gtexv8/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.lookup_table.txt.gz', sep="\t")

# get RSID and alleles for remaining SNPs
QTL.sig$rsid <- QTL.lookup$rs_id_dbSNP151_GRCh38p7[match(QTL.sig$variant_id, QTL.lookup$variant_id)] 
QTL.sig$ref <- QTL.lookup$ref[match(QTL.sig$variant_id, QTL.lookup$variant_id)]
QTL.sig$alt <- QTL.lookup$alt[match(QTL.sig$variant_id, QTL.lookup$variant_id)]
QTL.sig$locator <- paste(QTL.sig$gene_id, QTL.sig$rsid, QTL.sig$variant_id, sep="-")
QTL.sig$locator <- gsub("_b38", "", QTL.sig$locator)

# rm'd QTL so used this - count.variants <- 182641182
count.sig <- dim(QTL.sig)[1]
count.duplicated <- sum(duplicated(QTL.sig$rsid))

exposure_data <- format_data(dat = QTL.sig, type = "exposure",
                             snp_col = "locator",beta_col = "slope", 
                             se_col = "slope_se",eaf_col = "maf", 
                             effect_allele_col = "alt",other_allele_col = "ref",
                             pos_col = "locator",pval_col = "pval_nominal")

# for clumping, set "pos" as locator and SNP as "rsid" (hacky but package doesn't preserve custom columns)
exposure_data$pos <- exposure_data$SNP
exposure_data$SNP <- unlist(lapply(strsplit(exposure_data$pos, "-"), '[[' ,2))

# save processed data
write.csv(exposure_data, file = 'processed_GTExWholeBloodv8_eQTL.csv')
rm(QTL, QTL.sig, se)

# clump data 
exposure_data <- clump_data(dat = exposure_data, clump_kb = 10000, clump_r2 = 0.001)
write.csv(exposure_data, file = 'processed_GTExWholeBloodv8_eQTL.clumped.csv')

# get gene info
exposure_data$entrez <- toupper(unlist(lapply(strsplit(exposure_data$pos, "-"), '[[' ,1))) 
exposure_data$entrez <- unlist(lapply(strsplit(exposure_data$entrez, "[.]"), '[[', 1))

# snp and gene counts
count.snps <- length(exposure_data$SNP) 
count.snps.unique <- length(unique(exposure_data$SNP))
count.genes.unique <- length(unique(exposure_data$entrez))

rm(exposure_data)
rm(QTL.sig)
rm(count.duplicated, count.genes.unique, count.sig, count.snps, count.snps.unique, count.variants, 
   data_descr, individuals, individuals_descr, name,region, source)
rm(exposure_data)

write_metadata()
