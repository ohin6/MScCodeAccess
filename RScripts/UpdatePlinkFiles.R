###############
# Description #
###############
#* This script identified and corrects erroneous HLA genotyping from the .ped
#* file and converts back to plink --bfile format.
#* 
#* 
#* Steps:
#* 
#* 1. Install packages:
#* 2. In Plink: recode Plink files .bed/bim/fam to .ped/map format
#* 3. Import .ped and .map files as dataframe
#* 4. prepare .ped file
#* 5. Identify individuals with erroneous HLA genotyping
#* 6. Modify .ped file
#* 7. Re check for erroneous HLA genotyping
#* 8. Export new files
#* 9. In PLINK convert back to original PLINK --bfile format
#* 
#* Requirements:
#* Script requires original PLINK binary fileset HLA.subset.bed, HLA-subset.bim
#* and HLA-subset.fam saved in the same directory as the project.
#* 
#* Also an environment that allows PLINK to run through R, otherwise PLINK
#* operations will need to be performed outside this script.
#* 
#*  
## Author: Dr. Owen Williams
##
## Date Created: 05-12-2022
##
## Email: owen.williams8@nhs.net


#######################
# 1. Install packages #
#######################

require(trio)
require(tidyverse)

###################################################
# 2. In PLINK Recode bfiles into .ped/.map format #
###################################################
#* Using PLINK recode --bfiles into .ped/.map format
#* PLINK files need to be downloaded to current directory

system('plink --file Updated_HLA-subset --make-bed --out Updated_HLA-subset')

#* If PLINK environment is not set up the above code wont work and will need to 
#* be run in PLINK separately 
 
###################
# 3. Import files #
###################

# Import .ped file
#* Original .ped file with HLA genotyping
SNP_ped =  read.pedfile('../HLA_logitudinal_study/HLA_Plink_files/hla-subset.ped')

# Import SNP .map file
#* To get SNP names
SNP_map = tibble(read.table(file = '../HLA_logitudinal_study/HLA_Plink_files/hla-subset.map', sep = '\t'))

########################
# 4. Prepare .ped file #
########################

#Get .ped file colnames
SNPnames = SNP_map$V2
colnames_pedFile = colnames(SNP_ped[1:6])
                          
# loop function to get both alleles
for (i in SNPnames){
  x = paste0(i,"a")
  y = paste0(i,"b")
  colnames_pedFile = append(colnames_pedFile, x)
  colnames_pedFile = append(colnames_pedFile, y)
}

# apply column names
colnames(SNP_ped) = colnames_pedFile

#########################################################
# 5. Identify individuals with erroneous HLA genotyping #
#########################################################
#* The expected present/absent alleles should be 32/816. Any individuals outside
#* this range is considered erroneous.

# identif number of present/absent allele per individual
Present = rowSums(SNP_ped[7:854] == 'P')
Absent = rowSums(SNP_ped[7:854] == 'A')

# Create table showing individual wit erroneous HLA gebotyping
erroneous_individuals = SNP_ped %>%
  mutate(Number_Present_allele = as.numeric(Present)) %>%
  mutate(Number_Absent_allele = as.numeric(Absent)) %>%
  mutate(id = row_number()) %>%
  select(id, famid, Number_Present_allele, Number_Absent_allele) %>%
  filter(Number_Present_allele != 32 | Number_Absent_allele != 816)

#  Show table
erroneous_individuals
# 

#######################
# 6. Modify .ped file #
####################### 

# Create a row of 'Absent alleles'
vec <- setNames(rep("A", 848), colnames_pedFile[7:854])

# Create updated HLA genotyping rows
Update_id53 = bind_rows(vec)[1, ] %>%
  mutate(across(c("HLA_A_0101a", "HLA_A_0102a", "HLA_C_0501a", "HLA_C_1202a", "HLA_B_4402a", "HLA_B_5201a",
                  "HLA_DRB1_1301a","HLA_DRB1_1502a","HLA_DQA1_0103a","HLA_DQA1_0103b", "HLA_DQB1_0601a",
                  "HLA_DQB1_0603a","HLA_DPA1_0103a","HLA_DPA1_0201a","HLA_DPB1_0401a","HLA_DPB1_1701a",
                  "HLA_A_01a", "HLA_A_01b", "HLA_C_05a", "HLA_C_12a", "HLA_B_44a", "HLA_B_52a",
                  "HLA_DRB1_13a","HLA_DRB1_15a","HLA_DQA1_01a","HLA_DQA1_01b", "HLA_DQB1_06a",
                  "HLA_DQB1_06b","HLA_DPA1_01a","HLA_DPA1_02a","HLA_DPB1_04a","HLA_DPB1_17a"),
                ~replace(.x, .x == 'A', 'P')))
Update_id54 = bind_rows(vec)[1, ] %>%
  mutate(across(c("HLA_A_0201a", "HLA_A_0301a", "HLA_C_0401a", "HLA_C_1402a", "HLA_B_1516a","HLA_B_3501a",
                   "HLA_DRB1_0101a","HLA_DRB1_0102b","HLA_DQA1_0101a","HLA_DQA1_0101b", "HLA_DQB1_0501a",
                  "HLA_DQB1_0501b","HLA_DPA1_0103a","HLA_DPA1_0103b","HLA_DPB1_0401a", "HLA_DPB1_0402a",
                  "HLA_A_02a", "HLA_A_03a", "HLA_C_04a", "HLA_C_14a", "HLA_B_15a","HLA_B_35a",
                  "HLA_DRB1_01a","HLA_DRB1_01b","HLA_DQA1_01a","HLA_DQA1_01b", "HLA_DQB1_05a",
                  "HLA_DQB1_05b","HLA_DPA1_01a","HLA_DPA1_01b","HLA_DPB1_04a", "HLA_DPB1_04b"),
                ~replace(.x, .x == 'A', 'P')))
Update_id61 = bind_rows(vec)[1, ] %>%
  mutate(across(c('HLA_A_0201a', 'HLA_A_2402a', 'HLA_C_0304a', 'HLA_C_0305a', 'HLA_B_1401a', 'HLA_B_4001a',
                    'HLA_DRB1_0701a', 'HLA_DRB1_0701b', 'HLA_DQA1_0201a', 'HLA_DQA1_0201b', 'HLA_DQB1_0202a',
                    'HLA_DQB1_0202b', 'HLA_DPA1_0103a', 'HLA_DPA1_0103b', 'HLA_DPB1_0201a', 'HLA_DPB1_0201b',
                  'HLA_A_02a', 'HLA_A_24a', 'HLA_C_03a', 'HLA_C_03b', 'HLA_B_14a', 'HLA_B_40a',
                  'HLA_DRB1_07a', 'HLA_DRB1_07b', 'HLA_DQA1_02a', 'HLA_DQA1_02b', 'HLA_DQB1_02a',
                  'HLA_DQB1_02b', 'HLA_DPA1_01a', 'HLA_DPA1_01b', 'HLA_DPB1_02a', 'HLA_DPB1_02b'),
                ~replace(.x, .x == 'A', 'P')))

Update_id217 = bind_rows(vec)[1, ] %>%
  mutate(across(c('HLA_A_0201a', 'HLA_A_2601a', 'HLA_C_0702a', 'HLA_C_0702b', 'HLA_B_0702a', 'HLA_B_4102a',
                  'HLA_DRB1_0101a', 'HLA_DRB1_1303a', 'HLA_DQA1_0101a', 'HLA_DQA1_0501a', 'HLA_DQB1_0201a',
                  'HLA_DQB1_0501a', 'HLA_DPA1_0103a', 'HLA_DPA1_0201a', 'HLA_DPB1_0101a', 'HLA_DPB1_0201a',
                  'HLA_A_02a', 'HLA_A_26a', 'HLA_C_07a', 'HLA_C_07b', 'HLA_B_07a', 'HLA_B_41a',
                  'HLA_DRB1_01a', 'HLA_DRB1_13a', 'HLA_DQA1_01a', 'HLA_DQA1_05a', 'HLA_DQB1_02a',
                  'HLA_DQB1_05a', 'HLA_DPA1_01a', 'HLA_DPA1_02b', 'HLA_DPB1_01a', 'HLA_DPB1_02a'),
                ~replace(.x, .x == 'A', 'P')))
Update_id399 = bind_rows(vec)[1, ] %>%
  mutate(across(c('HLA_A_0101a', 'HLA_A_1101a', 'HLA_C_0501a', 'HLA_C_1202a', 'HLA_B_4402a', 'HLA_B_5201a',
                  'HLA_DRB1_1301a', 'HLA_DRB1_1502a', 'HLA_DQA1_0103a', 'HLA_DQA1_0103b', 'HLA_DQB1_0601a',
                  'HLA_DQB1_0603a', 'HLA_DPA1_0103a', 'HLA_DPA1_0201a', 'HLA_DPB1_0401a', 'HLA_DPB1_1701a',
                  'HLA_A_01a', 'HLA_A_11a', 'HLA_C_05a', 'HLA_C_12a', 'HLA_B_44a', 'HLA_B_52a',
                  'HLA_DRB1_13a', 'HLA_DRB1_15a', 'HLA_DQA1_01a', 'HLA_DQA1_01b', 'HLA_DQB1_06a',
                  'HLA_DQB1_06b', 'HLA_DPA1_01a', 'HLA_DPA1_02a', 'HLA_DPB1_04a', 'HLA_DPB1_17a'),
                ~replace(.x, .x == 'A', 'P')))
Update_id696 = bind_rows(vec)[1, ] %>%
  mutate(across(c('HLA_A_0201a','HLA_A_1101a','HLA_C_0701a','HLA_C_0702a','HLA_B_1801a','HLA_B_1801b',
                  'HLA_DRB1_0401a','HLA_DRB1_0404a','HLA_DQA1_0301a','HLA_DQA1_0301b','HLA_DQB1_0302a',
                  'HLA_DQB1_0302b', 'HLA_DPA1_0103a', 'HLA_DPA1_0201a', 'HLA_DPB1_1101a', 'HLA_DPB1_2301a',
                  'HLA_A_02a','HLA_A_11a','HLA_C_07a','HLA_C_07b','HLA_B_18a','HLA_B_18b',
                  'HLA_DRB1_04a','HLA_DRB1_04b','HLA_DQA1_03a','HLA_DQA1_03b','HLA_DQB1_03a',
                  'HLA_DQB1_03b', 'HLA_DPA1_01a', 'HLA_DPA1_02a', 'HLA_DPB1_11a', 'HLA_DPB1_23a'),
                ~replace(.x, .x == 'A', 'P')))
Update_id714 = bind_rows(vec)[1, ] %>%
  mutate(across(c('HLA_A_0101a', 'HLA_A_1101a', 'HLA_C_0202a', 'HLA_C_0701a', 'HLA_B_0702a', 'HLA_B_0801a',
                  'HLA_DRB1_0301a', 'HLA_DRB1_0701a', 'HLA_DQA1_0201a', 'HLA_DQA1_0501a', 'HLA_DQB1_0201a',
                  'HLA_DQB1_0303a', 'HLA_DPA1_0103a', 'HLA_DPA1_0201a', 'HLA_DPB1_0101a', 'HLA_DPB1_0401a',
                  'HLA_A_01a', 'HLA_A_11a', 'HLA_C_02a', 'HLA_C_07a', 'HLA_B_07a', 'HLA_B_08a',
                  'HLA_DRB1_03a', 'HLA_DRB1_07a', 'HLA_DQA1_02a', 'HLA_DQA1_05a', 'HLA_DQB1_02a',
                  'HLA_DQB1_03a', 'HLA_DPA1_01a', 'HLA_DPA1_02a', 'HLA_DPB1_01a', 'HLA_DPB1_04a'),
                ~replace(.x, .x == 'A', 'P')))
Update_id849 = bind_rows(vec)[1, ] %>%
  mutate(across(c('HLA_A_0201a', 'HLA_A_2902a', 'HLA_C_0602a', 'HLA_C_1601a', 'HLA_B_4404a', 'HLA_B_5701a',
                  'HLA_DRB1_0701a', 'HLA_DRB1_0701b', 'HLA_DQA1_0201a', 'HLA_DQA1_0201b', 'HLA_DQB1_0202a',
                  'HLA_DQB1_0303a', 'HLA_DPA1_0103a', 'HLA_DPA1_0103b', 'HLA_DPB1_0301a', 'HLA_DPB1_0401a',
                  'HLA_A_02a', 'HLA_A_29a', 'HLA_C_06a', 'HLA_C_16a', 'HLA_B_44a', 'HLA_B_57a',
                  'HLA_DRB1_07a', 'HLA_DRB1_07b', 'HLA_DQA1_02a', 'HLA_DQA1_02b', 'HLA_DQB1_02a',
                  'HLA_DQB1_03a', 'HLA_DPA1_01a', 'HLA_DPA1_01b', 'HLA_DPB1_03a', 'HLA_DPB1_04a'),
                ~replace(.x, .x == 'A', 'P')))
Update_id970 = bind_rows(vec)[1, ] %>%
  mutate(across(c('HLA_A_0301a', 'HLA_A_6802a', 'HLA_C_0401a', 'HLA_C_0702a', 'HLA_B_0702b', 'HLA_B_5301a',
                  'HLA_DRB1_0401a', 'HLA_DRB1_1302a', 'HLA_DQA1_0102a', 'HLA_DQA1_0301a', 'HLA_DQB1_0604a',
                  'HLA_DQB1_0604b', 'HLA_DPA1_0103a', 'HLA_DPA1_0103b', 'HLA_DPB1_0401a', 'HLA_DPB1_0401b',
                  'HLA_A_03a', 'HLA_A_68a', 'HLA_C_04a', 'HLA_C_07a', 'HLA_B_07b', 'HLA_B_53a',
                  'HLA_DRB1_04a', 'HLA_DRB1_13a', 'HLA_DQA1_01a', 'HLA_DQA1_03a', 'HLA_DQB1_06a',
                  'HLA_DQB1_06b', 'HLA_DPA1_01a', 'HLA_DPA1_01b', 'HLA_DPB1_04a', 'HLA_DPB1_04b'),
                ~replace(.x, .x == 'A', 'P')))
Update_id1275 = bind_rows(vec)[1, ] %>%
  mutate(across(c('HLA_A_0201a', 'HLA_A_2402a', 'HLA_C_0202a', 'HLA_C_0303a', 'HLA_B_1507a', 'HLA_B_2705a',
                    'HLA_DRB1_0404a', 'HLA_DRB1_0404b', 'HLA_DQA1_0101a', 'HLA_DQA1_0301a', 'HLA_DQB1_0302a',
                    'HLA_DQB1_0501a', 'HLA_DPA1_0103a', 'HLA_DPA1_0201a', 'HLA_DPB1_0601a', 'HLA_DPB1_0901a',
                    'HLA_A_02a', 'HLA_A_24a', 'HLA_C_02a', 'HLA_C_03a', 'HLA_B_15a', 'HLA_B_27a',
                    'HLA_DRB1_04a', 'HLA_DRB1_04b', 'HLA_DQA1_01a', 'HLA_DQA1_03a', 'HLA_DQB1_03a',
                    'HLA_DQB1_05a', 'HLA_DPA1_01a', 'HLA_DPA1_02a', 'HLA_DPB1_06a', 'HLA_DPB1_09a'),
                ~replace(.x, .x == 'A', 'P')))
Update_id1370 = bind_rows(vec)[1, ] %>%
  mutate(across(c('HLA_A_0301a', 'HLA_A_0301b', 'HLA_C_0102a', 'HLA_C_0401a', 'HLA_B_3501a', 'HLA_B_5101a',
                  'HLA_DRB1_0101a', 'HLA_DRB1_0101b', 'HLA_DQA1_0101a', 'HLA_DQA1_0101b', 'HLA_DQB1_0501a',
                  'HLA_DQB1_0501b', 'HLA_DPA1_0103a', 'HLA_DPA1_0103b', 'HLA_DPB1_0301a', 'HLA_DPB1_0402a',
                  'HLA_A_03a', 'HLA_A_03b', 'HLA_C_01a', 'HLA_C_04a', 'HLA_B_35a', 'HLA_B_51a',
                  'HLA_DRB1_01a', 'HLA_DRB1_01b', 'HLA_DQA1_01a', 'HLA_DQA1_01b', 'HLA_DQB1_05a',
                  'HLA_DQB1_05b', 'HLA_DPA1_01a', 'HLA_DPA1_01b', 'HLA_DPB1_03a', 'HLA_DPB1_04a'),
                ~replace(.x, .x == 'A', 'P')))
Update_id1423 = bind_rows(vec)[1, ] %>%
  mutate(across(c('HLA_A_0101a', 'HLA_A_0101b', 'HLA_C_0501a', 'HLA_C_0701a', 'HLA_B_0801a', 'HLA_B_4402a',
                  'HLA_DRB1_0301a', 'HLA_DRB1_0402a', 'HLA_DQA1_0301a', 'HLA_DQA1_0501a', 'HLA_DQB1_0201a',
                  'HLA_DQB1_0201b', 'HLA_DPA1_0103a', 'HLA_DPA1_0103b', 'HLA_DPB1_0201a', 'HLA_DPB1_0301a',
                  'HLA_A_01a', 'HLA_A_01b', 'HLA_C_05a', 'HLA_C_07a', 'HLA_B_08a', 'HLA_B_44a',
                  'HLA_DRB1_03a', 'HLA_DRB1_04a', 'HLA_DQA1_03a', 'HLA_DQA1_05a', 'HLA_DQB1_02a',
                  'HLA_DQB1_02b', 'HLA_DPA1_01a', 'HLA_DPA1_01b', 'HLA_DPB1_02a', 'HLA_DPB1_03a'),
                ~replace(.x, .x == 'A', 'P')))


# Change rows in .ped file
SNP_ped[53,][7:854] = Update_id53
SNP_ped[54,][7:854] = Update_id54
SNP_ped[61,][7:854] = Update_id61
SNP_ped[217,][7:854] = Update_id217
SNP_ped[399,][7:854] = Update_id399
SNP_ped[696,][7:854] = Update_id696
SNP_ped[714,][7:854] = Update_id714
SNP_ped[849,][7:854] = Update_id849
SNP_ped[970,][7:854] = Update_id970
SNP_ped[1275,][7:854] = Update_id1275
SNP_ped[1370,][7:854] = Update_id1370
SNP_ped[1423,][7:854] = Update_id1423

############################################
# 7. Re check for erroneous HLA genotyping #
############################################
#* Confirm all individuals have normal range of present and absent alleles
#* (32/816)

# Count number of present and absent
Present = rowSums(SNP_ped[7:854] == 'P')
Absent = rowSums(SNP_ped[7:854] == 'A')

# Create table showing individual wit erroneous HLA gebotyping
#* If table empty then all individual have normal range of alleles
ReCheckerroneous_individuals = SNP_ped %>%
  mutate(Number_Present_allele = as.numeric(Present)) %>%
  mutate(Number_Absent_allele = as.numeric(Absent)) %>%
  mutate(id = row_number()) %>%
  select(id, famid, Number_Present_allele, Number_Absent_allele) %>%
  filter(Number_Present_allele != 32 | Number_Absent_allele != 816)

###########################
# 8. Export new .ped file #
###########################

# Export corrected .ped file
write.table(SNP_ped, '../HLA_logitudinal_study/HLA_Plink_files/updated_HLA_alleleTyping/Updated_HLA-subset.ped',
            sep = " ", row.names = FALSE, col.names = FALSE, quote=FALSE)

# Export .map file with new name
write.table(SNP_map, '../HLA_logitudinal_study/HLA_Plink_files/updated_HLA_alleleTyping/Updated_HLA-subset.map',
            sep = " ", row.names = FALSE, col.names = FALSE, quote=FALSE)


############################################################
# 9. In Plink Convert .ped/.map file back to binary format #
############################################################
#* Using Plink 1.9 recode .ped and .map files back to binary format. This
#* requires running R in a preset environment. Otherwise Run below scripts in
#* PLINK seperately

system('plink --file Updated_HLA-subset --make-bed --out Updated_HLA-subset')