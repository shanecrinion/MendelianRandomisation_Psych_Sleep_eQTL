## Script to analysis eQTLs by the individual brain regions. 


# Libraries
library(dplyr)
library(TwoSampleMR)
library(rtracklayer)

# 1. Cortex

# Import and harmonise exposure and outcome data
exposure_data <- read.csv('processed_GTExCortex_eQTL.clumped.csv')
exposure_data$SNP <- unlist(lapply(strsplit(exposure_data$pos, "-"), '[[' ,2))
outcome_data <- extract_outcome_data(snps = exposure_data$SNP, outcomes = c('ieu-a-1183', 'ieu-a-1185', "ieu-b-41", "finn-b-F5_INSOMNIA", "ieu-b-102", "ukb-b-4956", "finn-b-F5_SCHZPHR"))
harm_data <- harmonise_data(exposure_data, outcome_data)
harm_data$SNP <- harm_data$pos.y

# Perform MR (wald ratio)
mr_wr <- mr_singlesnp(dat = harm_data)
mr_wald_ratio(singlesnp$beta.exposure, singlesnp$beta.outcome, 
              singlesnp$se.exposure, singlesnp$se.outcome)

# Extract significant results
mr_wr.sig <- subset(mr_wr, subset = p < 0.05) %>% arrange(p)
mr_wr.sig$entrez<- toupper(unlist(lapply(strsplit(mr_wr.sig$SNP, "-"), '[[' ,1)))
mr_wr.sig$entrez_base <- toupper(unlist(lapply(strsplit(mr_wr.sig$entrez, "[.]"), '[[' ,1)))
mr_wr.sig$rsid  <-unlist(lapply(strsplit(mr_wr.sig$SNP, "-"), '[[' ,2))

# Import gene names for significant genes
my_gtf <- import('gencode.v19.genes.v7.patched_contigs.gtf')
my_gtf.df<- my_gtf 
mr_wr.sig$gene_name <-  my_gtf.df$gene_name[match(mr_wr.sig$entrez, my_gtf.df$gene_id)]
dim(subset(mr_wr, p<(0.05/dim(mr_wr)[1])))

# Write results
write.csv(subset(mr_wr.sig, p<(0.05/4778)) %>% arrange(p),
          file = 'mr_wr_GTExCortex.alloutcomes.Bonf.csv')
write.csv(mr_wr.sig %>% arrange(p),file = 'mr_wr_GTExCortex.alloutcomes.p005.csv')
write.csv(mr_wr.sig %>% arrange(outcome,p),file = 'mr_wr_GTExCortex.byoutcome.alloutcomes.csv')

rm(exposure_data, harm_data, mr_wr, mr_wr.sig, outcome_data)

# 2. Blood

# Import and harmonise exposure and outcome data
exposure_data <- read.csv('processed_GTExWholeBlood_eQTL.clumped.csv')
exposure_data$SNP <- unlist(lapply(strsplit(exposure_data$pos, "-"), '[[' ,2))
outcome_data <- extract_outcome_data(snps = exposure_data$SNP, outcomes = c('ieu-a-1183', 'ieu-a-1185', "ieu-b-41", "finn-b-F5_INSOMNIA", "ieu-b-102", "ukb-b-4956", "finn-b-F5_SCHZPHR"))
harm_data <- harmonise_data(exposure_data, outcome_data)
harm_data$SNP <- harm_data$pos.y

# Perform MR (Wald ratio)
mr_wr <- mr_singlesnp(dat = harm_data)
mr_wr.sig <- subset(mr_wr, subset = p < 0.05) %>% arrange(p)

# Assign gene names
mr_wr.sig$entrez<- toupper(unlist(lapply(strsplit(mr_wr.sig$SNP, "-"), '[[' ,1)))
mr_wr.sig$entrez_base <- toupper(unlist(lapply(strsplit(mr_wr.sig$entrez, "[.]"), '[[' ,1)))
mr_wr.sig$rsid  <-unlist(lapply(strsplit(mr_wr.sig$SNP, "-"), '[[' ,2))
mr_wr.sig$gene_name <-  my_gtf.df$gene_name[match(mr_wr.sig$entrez, my_gtf.df$gene_id)]

# Write results
write.csv(subset(mr_wr.sig, p<0.05/dim(mr_wr)[1]) %>% arrange(p),
          file = 'mr_wr_GTExBlood.alloutcomes.Bonf.csv')
write.csv(mr_wr.sig %>% arrange(p),file = 'mr_wr_GTExWholeBlood.alloutcomes.p005.csv')
write.csv(mr_wr.sig %>% arrange(outcome,p),file = 'mr_wr_GTExWholeBlood.byoutcome.alloutcomes.csv')

# 3. Hypothalamus

# Import and harmonise exposure and outcome data
exposure_data <- read.csv('processed_GTExHypothalamus_eQTL.clumped.csv')
exposure_data$SNP <- unlist(lapply(strsplit(exposure_data$pos, "-"), '[[' ,2))
outcome_data <- extract_outcome_data(snps = exposure_data$SNP, 
                                     outcomes = c('ieu-a-1183', 'ieu-a-1185', "ieu-b-41", 
                                     "finn-b-F5_INSOMNIA", "ieu-b-102", "ukb-b-4956",
                                      "finn-b-F5_SCHZPHR"))
harm_data <- harmonise_data(exposure_data, outcome_data)
harm_data$SNP <- harm_data$pos.y

# Perform MR
mr_wr <- mr_singlesnp(dat = harm_data)
mr_wr.sig <- subset(mr_wr, subset = p < 0.05) %>% arrange(p)

# Assign gene names
mr_wr.sig$entrez<- toupper(unlist(lapply(strsplit(mr_wr.sig$SNP, "-"), '[[' ,1)))
mr_wr.sig$entrez_base <- toupper(unlist(lapply(strsplit(mr_wr.sig$entrez, "[.]"), '[[' ,1)))
mr_wr.sig$rsid  <-unlist(lapply(strsplit(mr_wr.sig$SNP, "-"), '[[' ,2))
mr_wr.sig$gene_name <-  my_gtf.df$gene_name[match(mr_wr.sig$entrez, my_gtf.df$gene_id)]

# Write results
write.csv(subset(mr_wr.sig, p<0.05/dim(mr_wr)[1]) %>% arrange(p),
          file = 'mr_wr_GTExHypothalamus.alloutcomes.Bonf.csv')
write.csv(mr_wr.sig %>% arrange(p),file = 'mr_wr_GTExHypothalamus.alloutcomes.p005.csv')
write.csv(mr_wr.sig %>% arrange(outcome,p),file = 'mr_wr_GTExHypothalamus.byoutcome.alloutcomes.csv')

