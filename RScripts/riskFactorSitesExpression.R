###############
# Description #
###############
#* This script examines how the identified risk factor sites effect gene 
#* expression in microglia tissue.
#* 
#* eQTL data examines gene expresssion in:
#*  1. HLA-DRB1 (ENSG00000196126)
#*  2. HLA-DQA1 (ENSG00000196735)
#*  3. HLA-DRA (ENSG00000204287)
#* 
#* Steps
#* 1. Set up link with online eQTL catalog
#* 2. get eQTL data
#* 3. Convert to GRCH37 build
#* 4. Get proportions of increased and decreased gene expression
#* 5. Logistic regression Model to see whether region risk factor regions reduce risk expression levels
#* 
#* 
## Author: Dr. Owen Williams
##
## Date Created: 02-01-2023
##
## Email: owen.williams8@nhs.net


####################
# Install packages #
####################

require(tidyverse)
library("readr")
library("seqminer")
require(rms)

###########################################
# 1. Set up link with online eQTL catalog #
###########################################

# Set link to online eQTL catalog
tabix_paths = read.delim("https://raw.githubusercontent.com/eQTL-Catalogue/eQTL-Catalogue-resources/master/tabix/tabix_ftp_paths.tsv", sep = "\t", header = TRUE, stringsAsFactors = FALSE) %>% dplyr::as_tibble()

# Get list of eQTL data from GTex study 
imported_tabix_paths = read.delim("https://raw.githubusercontent.com/eQTL-Catalogue/eQTL-Catalogue-resources/master/tabix/tabix_ftp_paths_imported.tsv", sep = "\t", header = TRUE, stringsAsFactors = FALSE) %>% dplyr::as_tibble()

# Create function for extracting information from study using Tabix
import_eQTLCatalogue <- function(ftp_path, region, selected_gene_id, column_names, verbose = TRUE){
  
  if(verbose){
    print(ftp_path)
  }
  
  #Fetch summary statistics with seqminer
  fetch_table = seqminer::tabix.read.table(tabixFile = ftp_path, tabixRange = region, stringsAsFactors = FALSE) %>%
    dplyr::as_tibble()
  colnames(fetch_table) = column_names
  
  #Remove rsid duplicates and multi-allelic variant
  summary_stats = dplyr::filter(fetch_table, gene_id == selected_gene_id) %>%
    dplyr::select(-rsid) %>% 
    dplyr::distinct() %>% #rsid duplicates
    dplyr::mutate(id = paste(chromosome, position, sep = ":")) %>% 
    dplyr::group_by(id) %>% 
    dplyr::mutate(row_count = n()) %>% dplyr::ungroup() %>% 
    dplyr::filter(row_count == 1) #Multialllics
}

###################
# 2. Extract eQTL #
###################
#* Region of interest is the HLA QTL on chromosome 6. The GTEx uses the GRCh38 build
#* so this will need to be converted to GRCh37 build afterwards. The region of
#* interest in GRCh37 is Chr6:32300000-32750000 which when converted to GRCh38 
#* is Chr6:32332223-32782223

# Define region of interest in GRCh38 build
region = "6:32332223-32782223"


#identify studies in GTEx
Microglia = dplyr::filter(tabix_paths, study == "Young_2019")

# 3.1 Create function -----------------------------------------------------

#* Create function plot expression results HLA region in user selected
#* tissue type (see above to see available tissues)
HLA_expression_tissue = function(tissue, selectedGene){
  studies_df = dplyr::filter(tabix_paths, study == "Young_2019", tissue_label == tissue)
  #Extract column names from first file
  column_names = colnames(readr::read_tsv(studies_df$ftp_path, n_max = 1))
  #Import summary statistics - selected gene = 
  summary_stats1 = import_eQTLCatalogue(studies_df$ftp_path, region, selected_gene_id = selectedGene, column_names) 
  # plot
  plot1 = ggplot(summary_stats1, aes(x = position, y = -log(pvalue, 10))) +
    labs(title = tissue) +
    geom_point() +
    ylim(0,max(-log(summary_stats1$pvalue,10)))
  list = list(plot1, summary_stats1)
  return(list)
}

# Extract eQTL -----------------------------------------------------------
#* Get ensemble gene names
DRB1 = "ENSG00000196126"
DQA1 = "ENSG00000196735"
DRA = "ENSG00000204287"

variants = c(32431147, 32441199, 32544901, 32550331, 32550828, 32550860, 32560631, 32561331,
32572113, 32572302, 32572305, 32575513, 32576592, 32576688, 32577222, 32578772,
32583027, 32583357, 32584693, 32606941, 32606949, 32607141, 32607729, 32610825,
32611590, 32619144) 


#* Examine HLA-DRB1 expression in microglia (repeat for others genes in subsequent analysis)
microglia_DRB1 = HLA_expression_tissue("microglia", DRB1)
microglia_DQA1 = HLA_expression_tissue("microglia", DQA1)
microglia_DRA = HLA_expression_tissue("microglia", DRA)


##############################
# 3. Convert to GRCh37 build #
##############################

# Convert column names to standardised format
sumStatsDRB1 = microglia_DRB1[[2]] %>%
  rename(CHR = chromosome) %>%
  rename(SNP = variant) %>%
  rename(BP = position) %>%
  rename(P = pvalue)

# Convert to GRCh38 to GRCh37 using mungestats lift over
sumStatsDRB1 = MungeSumstats::liftover(sumstats_dt = sumStatsDRB1,
                                   ref_genome = "hg38",
                                   convert_ref_genome = "hg19")


# Convert column names to standardised format
sumStatsDQA1 = microglia_DQA1[[2]] %>%
  rename(CHR = chromosome) %>%
  rename(SNP = variant) %>%
  rename(BP = position) %>%
  rename(P = pvalue)

# Convert to GRCh38 to GRCh37 using mungestats lift over
sumStatsDQA1 = MungeSumstats::liftover(sumstats_dt = sumStatsDQA1,
                                       ref_genome = "hg38",
                                       convert_ref_genome = "hg19")


# Convert column names to standardised format
sumStatsDRA = microglia_DRA[[2]] %>%
  rename(CHR = chromosome) %>%
  rename(SNP = variant) %>%
  rename(BP = position) %>%
  rename(P = pvalue)

# Convert to GRCh38 to GRCh37 using mungestats lift over
sumStatsDRA = MungeSumstats::liftover(sumstats_dt = sumStatsDRA,
                                       ref_genome = "hg38",
                                       convert_ref_genome = "hg19")

microglia_DRB1[[1]]
microglia_DQA1[[1]]
microglia_DRA[[1]]

#################################################################
# 4. Get proportions of increased and decreased gene expression #
#################################################################

# Risk factor site 1 [+/-500bp]
#* DRB1
sumStatsDRB1_RF1 =  sumStatsDRB1 %>%
  filter(between(BP, 32440699,32441699))
down = sum(sumStatsDRB1_RF1$beta < 0)
up = sum(sumStatsDRB1_RF1$beta >= 0)
proportion = (down/(down+up))*100
proportion

#*DQA1
sumStatsDQA1_RF1 =  sumStatsDQA1 %>%
  filter(between(BP, 32440699,32441699))
down = sum(sumStatsDQA1_RF1$beta < 0)
up = sum(sumStatsDQA1_RF1$beta >= 0)
proportion = (down/(down+up))*100
proportion

#*DRA
sumStatsDRA_RF1 =  sumStatsDRA %>%
  filter(between(BP, 32440699,32441699))
down = sum(sumStatsDRA_RF1$beta < 0)
up = sum(sumStatsDRA_RF1$beta >= 0)
proportion = (down/(down+up))*100
proportion


# Risk factor site 2
#*DRB1
sumStatsDRB1_RF2 = sumStatsDRB1 %>%
  filter(between(BP, 32571113,32573305))
down = sum(sumStatsDRB1_RF2$beta < 0)
up = sum(sumStatsDRB1_RF2$beta >= 0)
proportion = (down/(down+up))*100
proportion

#*DQA1
sumStatsDQA1_RF2 =  sumStatsDQA1 %>%
  filter(between(BP, 32571113,32573305))
down = sum(sumStatsDQA1_RF2$beta < 0)
up = sum(sumStatsDQA1_RF2$beta >= 0)
proportion = (down/(down+up))*100
proportion

#*DRA
sumStatsDRA_RF2 =  sumStatsDRA %>%
  filter(between(BP, 32571113,32573305))
down = sum(sumStatsDRA_RF2$beta < 0)
up = sum(sumStatsDRA_RF2$beta >= 0)
proportion = (down/(down+up))*100
proportion

# Risk factor site 3
#*DRB1
sumStatsDRB1_RF3 = sumStatsDRB1 %>%
  filter(between(BP, 32575513,32578772))
down = sum(sumStatsDRB1_RF3$beta < 0)
up = sum(sumStatsDRB1_RF3$beta >= 0)
proportion = (down/(down+up))*100
proportion

#*DQA1
sumStatsDQA1_RF3 =  sumStatsDQA1 %>%
  filter(between(BP, 32575513,32578772))
down = sum(sumStatsDQA1_RF3$beta < 0)
up = sum(sumStatsDQA1_RF3$beta >= 0)
proportion = (down/(down+up))*100
proportion

#*DRA
sumStatsDRA_RF3 =  sumStatsDRA %>%
  filter(between(BP, 32575513,32578772))
down = sum(sumStatsDRA_RF3$beta < 0)
up = sum(sumStatsDRA_RF3$beta >= 0)
proportion = (down/(down+up))*100
proportion


# All risk factor sites
#*DRB1
sumStatsDRB1_RF3 = sumStatsDRB1 %>%
  filter(between(BP, 32575513,32578772)|between(BP, 32571113,32573305)|between(BP, 32440699,32441699))
down = sum(sumStatsDRB1_RF3$beta < 0)
up = sum(sumStatsDRB1_RF3$beta >= 0)
proportion = (down/(down+up))*100
proportion

#*DQA1
sumStatsDQA1_RF3 =  sumStatsDQA1 %>%
  filter(between(BP, 32575513,32578772)|between(BP, 32571113,32573305)|between(BP, 32440699,32441699))
down = sum(sumStatsDQA1_RF3$beta < 0)
up = sum(sumStatsDQA1_RF3$beta >= 0)
proportion = (down/(down+up))*100
proportion

#*DRA
sumStatsDRA_RF3 =  sumStatsDRA %>%
  filter(between(BP, 32575513,32578772)|between(BP, 32571113,32573305)|between(BP, 32440699,32441699))
down = sum(sumStatsDRA_RF3$beta < 0)
up = sum(sumStatsDRA_RF3$beta >= 0)
proportion = (down/(down+up))*100
proportion



# input risk factor regions -----------------------------------------------
#* DRB1
sumStatsDRB1 = sumStatsDRB1 %>%
  mutate(riskFactor1 = ifelse(between(BP, 32440699,32441699), 1, 0)) %>%
  mutate(riskFactor2 = ifelse(between(BP, 32571113,32573305), 1, 0)) %>%
  mutate(riskFactor3 = ifelse(between(BP, 32575513,32578772), 1, 0)) %>%
  mutate(riskFactor = ifelse(between(BP, 32575513,32578772)|between(BP, 32571113,32573305)|between(BP, 32440699,32441699), 1, 0))

#* DQA1
sumStatsDQA1 = sumStatsDQA1 %>%
  mutate(riskFactor1 = ifelse(between(BP, 32440699,32441699), 1, 0)) %>%
  mutate(riskFactor2 = ifelse(between(BP, 32571113,32573305), 1, 0)) %>%
  mutate(riskFactor3 = ifelse(between(BP, 32575513,32578772), 1, 0)) %>%
  mutate(riskFactor = ifelse(between(BP, 32575513,32578772)|between(BP, 32571113,32573305)|between(BP, 32440699,32441699), 1, 0))

#* DRA
sumStatsDRA = sumStatsDRA %>%
  mutate(riskFactor1 = ifelse(between(BP, 32440699,32441699), 1, 0)) %>%
  mutate(riskFactor2 = ifelse(between(BP, 32571113,32573305), 1, 0)) %>%
  mutate(riskFactor3 = ifelse(between(BP, 32575513,32578772), 1, 0)) %>%
  mutate(riskFactor = ifelse(between(BP, 32575513,32578772)|between(BP, 32571113,32573305)|between(BP, 32440699,32441699), 1, 0))



# Export files ------------------------------------------------------------

write_csv(sumStatsDRB1, '../../write up/Supplementary/eQTLDRB1Microglia.csv')
write_csv(sumStatsDQA1, '../../write up/Supplementary/eQTLDQA1Microglia.csv')
write_csv(sumStatsDRA, '../../write up/Supplementary/eQTLDRAMicroglia.csv')



########################################################################################################
# 5. Logistic regression Model to see whether region risk factor regions reduce risk expression levels #
########################################################################################################

# Model for DRB1 ----------------------------------------------------------
# Get logistic regression for 
dd = datadist(sumStatsDRB1)
# describe distributions of variables to rms
options(datadist= 'dd' )

# Create model to see whether risk factor site reduce expression level
f_DRB1 = lrm(sumStatsDRB1$beta<0 ~ riskFactor1 + riskFactor2 + riskFactor3  x = TRUE, y = TRUE, data = sumStatsDRB1)
summary(f_DRB1)

# all risk factors combined
f_DRB1_com = lrm(sumStatsDRB1$beta<0 ~ riskFactor,  x = TRUE, y = TRUE, data = sumStatsDRB1)
summary(f_DRB1_com)

# Model for DQA1 ----------------------------------------------------------
# Get logistic regression for 
dd = datadist(sumStatsDQA1)
# describe distributions of variables to rms
options(datadist= 'dd' )

# Create model to see whether risk factor site reduce expression level
f_DQA1 = lrm(sumStatsDQA1$beta<0 ~ riskFactor1 + riskFactor2 + riskFactor3, x = TRUE, y = TRUE, data = sumStatsDQA1)
summary(f_DQA1)

# all risk factors combined
f_DQA1_com = lrm(beta<0 ~ riskFactor,  x = TRUE, y = TRUE, data = sumStatsDQA1)
summary(f_DQA1_com)

# Model for DRA ----------------------------------------------------------
# Get logistic regression for 
dd = datadist(sumStatsDRA)
# describe distributions of variables to rms
options(datadist= 'dd' )

# Create model to see whether risk factor site reduce expression level
f_DRA = lrm(sumStatsDRA$beta<0 ~ riskFactor1 + riskFactor2 + riskFactor3, x = TRUE, y = TRUE, data = sumStatsDRA)
summary(f_DRA)

# all risk factors combined
f_DRA_com = lrm(beta<0 ~ riskFactor,  x = TRUE, y = TRUE, data = sumStatsDRA)
summary(f_DRA_com)






