###############
# Description #
###############
#* The purpose of this script is to convert the SNP .ped file into a matrix which
#* can be used for later analysis, namely MiDAS
#* 
#* The .ped file comes from the cognitive decline longitudinal study provided by
#* Tony and comprises of different HLA haplotypes.
#* 
#* **NOTE** The .ped file was updated as the original had erroneous genotyping
#* 
#* 
## Author: Dr. Owen Williams
##
## Date Created: 07-12-2022
##
## Email: owen.williams8@nhs.net

#######################
# 1. Install packages #
#######################

require(tidyverse)

# set directory to current location
setwd(here::here())

###################
# 2. Import files #
###################

# SNP .ped file
#* gives SNPs per individual
SNP_ped =  tibble(read.table('../HLA_logitudinal_study/HLA_Plink_files/updated_HLA_alleleTyping/Updated_HLA-subset.ped'))

# SNP .map file
#* Gives SNP names
SNP_map = tibble(read.table(file = '../HLA_logitudinal_study/HLA_Plink_files/updated_HLA_alleleTyping/Updated_HLA-subset.map', sep = '\t'))

####################################
# 3. Convert SNPs to matrix format #
####################################

# Select patient ID and SNPs
SNPMatrix = SNP_ped %>%
  select(1, 7:854)

#* Because .ped file is in diploid format, create loop to get two copies of each
#* allele (a/b)
#* 

#Get row names from .map file
SNP_names = SNP_map$V2
col_SNPs= "ID"

# loop function to get both alleles
for (i in SNP_names){
  x = paste0(i,"a")
  y = paste0(i,"b")
  col_SNPs = append(col_SNPs, x)
  col_SNPs = append(col_SNPs, y)
}

# Add updated colnames
colnames(SNPMatrix) = col_SNPs

#####################################
# 4. Convert matrix to MiDAs format #
#####################################


# 4.1 Transpose dataframe -----------------------------------------------------

# transpose table SNP matrix table
Trans_Genotypes = as.data.frame(t(SNPMatrix))
#* Convert patient ID as colnames
colnames(Trans_Genotypes) = SNPMatrix$ID
Trans_Genotypes = Trans_Genotypes[-1,] # remove redundant patient id row


# Filter rows containing P ------------------------------------------------
#* Where P denote Present allele

# Get Patient Ids
patientIds = SNPMatrix$ID

# create empty tibble which contains different HLA genes

MiDASGeno = tibble("HLA_A1", "HLA_A2","HLA_A3","HLA_A4", "HLA_C1", "HLA_C2", 
                   "HLA_C3","HLA_C4", "HLA_B1", "HLA_B2","HLA_B3","HLA_B4",
                   "HLA_DRB1", "HLA_DRB12","HLA_DRB13","HLA_DRB14","HLA_DQA11",
                   "HLA_DQA12","HLA_DQA13","HLA_DQA14", "HLA_DQB11", "HLA_DQB12",
                   "HLA_DQB13","HLA_DQB14", "HLA_DPA11", "HLA_DPA12","HLA_DPA13",
                   "HLA_DPA14","HLA_DPB11", "HLA_DPB12","HLA_DPB13","HLA_DPB14")

# Create empty vector
dfRowNames = vector()

# For loop through each column in Trans_genotype dataframe which contains 'p' for present allele
for (i in seq_along(patientIds)){
  df = Trans_Genotypes[i] %>%
    filter(str_detect(Trans_Genotypes[[i]], "P"))
  dfRowNames = row.names(df)
  MiDASGeno = rbind(MiDASGeno,dfRowNames)
}

# Remove first row from matrix as this is redundant (same as colnames)
MiDASGeno = MiDASGeno[-1,]

# Change row names as patient IDs
rownames(MiDASGeno) = patientIds


# 4.3 Tidy Data -----------------------------------------------------------

# remove last letter from every element
for (i in 1:ncol(MiDASGeno)){
  MiDASGeno[[i]] = gsub("a$","", MiDASGeno[[i]])
  MiDASGeno[[i]] = gsub("b$","", MiDASGeno[[i]])
}

# Transpose MiDASGeno
MiDASGeno = as.data.frame(t(MiDASGeno))

# Change colnames of MiDASGeno to show allele types
colnames(MiDASGeno) = c("A_1", "A_2", "C_1", "C_2", "B_1", "B_2", "DRB1_1", 
                         "DRB1_2","DQA1_1", "DQA1_2", "DQB1_1", "DQB1_2", 
                         "DPA1_1", "DPA1_2","DPB1_1", "DPB1_2")

# clear data
MiDASGeno = MiDASGeno2
rm(MiDASGeno2)


# 4.7 Tidy data 2 -------------------------------------------------------------

# Order columns
MiDASGeno = MiDASGeno %>%
  select(A_1,	A_2,	B_1,	B_2,	C_1,	C_2,	DPA1_1,	DPA1_2,	DPB1_1,	DPB1_2,
         DQA1_1,	DQA1_2,	DQB1_1,	DQB1_2,	DRB1_1,	DRB1_2)

# Adjust nomenclature to match MiDAS
MiDASGeno = MiDASGeno %>%
  mutate(across(everything(), str_replace_all, 'HLA_', '')) %>%
  mutate(across(everything(), str_replace_all, '_', '*'))


# Add a ':' between last two digits of each HLA allele
#* Create function to add character 
fun_insert = function(x, pos, insert){
  gsub(paste0("^(.{", pos,"})(.*)$"),
       paste0("\\1",insert,"\\2"),
       x)
}

# Add : to second from last charcter of string
MiDASGeno = MiDASGeno %>%
  mutate(across(everything(), ~ fun_insert(.x,str_length(.x)-2,':')))

# Add ID column to start of matrix
MiDASGeno = MiDASGeno %>%
  add_column(ID = SNPMatrix$ID, .before = 'A_1')

########################
# 5. Export data table #
########################

# Export as .txt file
write_tsv(MiDASGeno, "../HLA_logitudinal_study/HLA_Plink_files/updated_HLA_alleleTyping/MiDASGeno.txt")




