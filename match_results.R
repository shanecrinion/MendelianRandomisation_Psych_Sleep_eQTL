# Rscript

# Author: Shane Crinion
# Contact: shanecrinion@gmail.com
# Description: R script to combine all MR and coloc results and create results files for each tissue type

# Function to format all results from coloc per SNP
match_coloc_results <- function(results_files, data){
  for(i in results_files){
    info <- unlist(strsplit(i, '_'))
    #print(unlist(strsplit(info[1], '/')))
    snp_i <- unlist(strsplit(info[1], '/'))[2]; print(snp_i)
    entrez_i <- info[2] ; print(entrez_i)
    outcome_i <- info[5]; print(outcome_i)
    snp_result <- subset(read.csv(paste0('coloc/', i)), snp==snp_i) [['SNP.PP.H4']]
    sum_result <- read.csv(stringr::str_replace(paste0('coloc/', i), "results", "summary"), row.names = 1)[6,]
    print(str(data))
    data[data$rsid==snp_i & data$entrez==entrez_i  & data$outcome==outcome_i,][['snp_pp4']] <- snp_result
    data[data$rsid==snp_i & data$entrez==entrez_i & data$outcome==outcome_i,][['sum_pp4']] <- sum_result
  }
  message("need to manually edit the following rows:") # Handling duplicates
  results=list(coloc=data, editing_needed=distinct(data[match(duplicated(data$rsid), 
                                                              duplicated(data$entrez)),]))
  return(results)}


# 1. ---- GTEx Cortex
# Read in results from MR analysis (to be combined w coloc results)
googlesheets4::gs4_auth(email = "shanecrinion@gmail.com")
data <-read_sheet(
  "https://docs.google.com/spreadsheets/d/1FygOe-3mZynCf1TvKIu_8F8EIfbxgndMpKW_EqbBN7Y/edit?usp=sharing",sheet = "GTEx v8 (cortex)")

# Extract RSID and Entrez ID to match to SNPs
data$rsid <- unlist(lapply(strsplit(data$`SNP (entrez-rsid-pos)`, "-"), '[[', 2))
data$entrez <- toupper(unlist(lapply(strsplit(data$`SNP (entrez-rsid-pos)`, "-"), '[[', 1)))

# Create empty values for the coloc results
data$snp_pp4 <- 0
data$sum_pp4 <- 0

# Get coloc results and write combined results file. 
results_files <- list.files('~/Desktop/files/coloc/', pattern = 'results.csv') # make recursive
results_files <- results_files[grepl("^rs.*_cortex_", results_files)]
cortex_results <- match_coloc_results(results_files = results_files, data = data)

write.csv(cortex_results$coloc, 'coloc/results_GTExCortexv8_eQTL_GWAS.Bonf.withcoloc.csv')


# 2. ----  GTEx Hippocampus
# Read in results from MR analysis (to be combined w coloc results)
googlesheets4::gs4_auth(email = "shanecrinion@gmail.com")
data <-read_sheet(
  "https://docs.google.com/spreadsheets/d/1FygOe-3mZynCf1TvKIu_8F8EIfbxgndMpKW_EqbBN7Y/edit?usp=sharing",sheet = "GTEx v8 (hippocampus)")

# Extract RSID and Entrez ID to match to SNPs
data$rsid <- unlist(lapply(strsplit(data$`SNP (entrez-rsid-pos)`, "-"), '[[', 2))
data$entrez <- toupper(unlist(lapply(strsplit(data$`SNP (entrez-rsid-pos)`, "-"), '[[', 1)))

# Create empty values for the coloc results
data$snp_pp4 <- 0
data$sum_pp4 <- 0

# Get coloc results and write combined results file. 
results_files <- list.files('~/Desktop/files/coloc/', pattern = 'results.csv')
results_files <- results_files[grepl("_hippocampus_", results_files)]
hippocampus_results <- match_coloc_results(results_files = results_files, data = data)
write.csv(hippocampus_results$coloc, 'coloc/results_GTExHippocampusv8_eQTL_GWAS.Bonf.withcoloc.csv')

# 3. ----  GTEx Hypothalamus
googlesheets4::gs4_auth(email = "shanecrinion@gmail.com")
data <-read_sheet(
  "https://docs.google.com/spreadsheets/d/1FygOe-3mZynCf1TvKIu_8F8EIfbxgndMpKW_EqbBN7Y/edit?usp=sharing",sheet = "GTEx v8 (hypothalamus)")

# Extract RSID and Entrez ID to match to SNPs
data$rsid <- unlist(lapply(strsplit(data$`SNP (entrez-rsid-pos)`, "-"), '[[', 2))
data$entrez <- toupper(unlist(lapply(strsplit(data$`SNP (entrez-rsid-pos)`, "-"), '[[', 1)))

# Create empty values for the coloc results
data$snp_pp4 <- 0
data$sum_pp4 <- 0

# Get coloc results and write combined results file. 
results_files <- list.files('~/Desktop/files/coloc/', pattern = 'beta_change.results.csv',recursive = T)
results_files <- results_files[grepl("_hypothalamus_", results_files)]
hypothalamus_results <- match_coloc_results(results_files = results_files, data = data)
write.csv(hypothalamus_results$coloc, '~/Desktop/files/coloc/results_GTExHypothalamusv8_eQTL_GWAS.Bonf.flip.withcoloc.csv')


# 4. ----  GTEx - Whole Blood 
data <-read_sheet(
  "https://docs.google.com/spreadsheets/d/1FygOe-3mZynCf1TvKIu_8F8EIfbxgndMpKW_EqbBN7Y/edit?usp=sharing",
  sheet = "GTEx (whole blood)")

# Extract RSID and Entrez ID to match to SNPs
data$rsid <- unlist(lapply(strsplit(data$`SNP (entrez-rsid-pos)`, "-"), '[[', 2))
data$entrez <- toupper(unlist(lapply(strsplit(data$`SNP (entrez-rsid-pos)`, "-"), '[[', 1)))

# Create empty values for the coloc results
data$snp_pp4 <- 0
data$sum_pp4 <- 0

# Get coloc results and write combined results file
results_files <- list.files('~/Desktop/files', pattern = 'results.csv')
results_files <- results_files[grepl("^rs.*_blood_", results_files)]
blood_results <-  match_coloc_results(results_files = results_files, data = data)
write.csv(blood_results$coloc, '~/Desktop/files/coloc/results_GTEx_wholeblood_GWAS.Bonf.withcoloc.csv')

# .  ----  PsychENCODE (prefrontal cortex)
googlesheets4::gs4_auth(email = "shanecrinion@gmail.com")
data <-read_sheet(
  "https://docs.google.com/spreadsheets/d/1FygOe-3mZynCf1TvKIu_8F8EIfbxgndMpKW_EqbBN7Y/edit?usp=sharing",sheet = "PsychENCODE (prefrontal cortex)")

# Extract RSID and Entrez ID to match to SNPs
data$rsid <- unlist(lapply(strsplit(data$`SNP (entrez-rsid-pos)`, "-"), '[[', 2))
data$entrez <- toupper(unlist(lapply(strsplit(data$`SNP (entrez-rsid-pos)`, "-"), '[[', 1)))

# Create empty values for the coloc results
data$snp_pp4 <- 0
data$sum_pp4 <- 0

# Get coloc results and write combined results file. 
results_files <- list.files('coloc/', pattern = 'results.csv', recursive = T)
results_files <- results_files[grepl("prefrontal", results_files)]
prefrontalcortex_results <- match_coloc_results(results_files = results_files, data = data)
write.csv(prefrontalcortex_results$coloc, 'coloc/results_PSYCHENCODE_eQTL_GWAS.Bonf.withcoloc.csv')




