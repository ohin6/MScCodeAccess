###############
# Description #
###############
#* The purpose of this script is to tidy RAW phenotypic data from cohort study 
#* into a format that will facilitate later association analysis in PLINK
#* 
#* Steps:
#* 
#* 1. Import Phenotype data
#* 2. Integrate APOE genotype data
#* 3. Tidy Data
#* 4. Impute missing data
#* 5. Examine data
#* 
#* 
## Author: Dr. Owen Williams
##
## Date Created: 07-12-2022
##
## Email: owen.williams8@nhs.net


####################
# Install packages #
####################

library(tidyverse)
library(readxl)
library(stringr)
require(here)
require(rms)
require(ggpubr)

#####################
# 1. Import dataset #
#####################

# set current directory
setwd(here::here())

# Import dataset
Brain_pathology <- tibble(read_excel("../Raw_data/CohortStudy/PhenotypicData/2021_11_16_Brain pathology.xlsx"))

###################################
# 2. Integrate APOE genotype data #
###################################

# Import APOE data
APOE = read_excel("../Raw_data/CohortStudy/PhenotypicData/APOE_E4pres_geno.xlsx") %>%
  select(panelid, apoe_e4, apoe)

# Join with Brain phenotype data
Brain_pathology = Brain_pathology %>%
  left_join(APOE, c("ID" = "panelid"))

# Remove APOE file to save memory
rm(APOE)

################
# 3. Tidy data #
################

# 1.  Convert BRAAK stage from roman numerals to an ordinal predictor ----------

#Convert to numeric
Brain_pathology = Brain_pathology %>%
  separate (braakstage, into= c("braakstage1", "braakstage2")) %>%
  mutate(braakstage1 = ifelse(braakstage1== "VI",6,braakstage1)) %>%
  mutate(braakstage1 = ifelse(braakstage1== "IV",4,braakstage1)) %>%
  mutate(braakstage1 = ifelse(braakstage1== "V",5,braakstage1)) %>%
  mutate(braakstage1 = ifelse(braakstage1== "III",3,braakstage1)) %>%
  mutate(braakstage1 = ifelse(braakstage1== "II",2,braakstage1)) %>%
  mutate(braakstage1 = ifelse(braakstage1== "I",1,braakstage1)) %>%
  mutate(braakstage1 = ifelse(braakstage1== "N", NA, braakstage1)) %>%
  mutate(braakstage2 = ifelse(braakstage2== "VI",6,braakstage2)) %>%
  mutate(braakstage2 = ifelse(braakstage2=="IV",4, braakstage2)) %>%
  mutate(braakstage2 = ifelse(braakstage2== "V",5, braakstage2)) %>%
  mutate(braakstage2 = ifelse(braakstage2== "III",3, braakstage2)) %>%
  mutate(braakstage2 = ifelse(braakstage2== "II",2, braakstage2)) %>%
  mutate(braakstage2 = ifelse(braakstage2== "I",1, braakstage2)) %>%
  mutate(braakstage2 = ifelse(braakstage2== "A", NA, braakstage2)) %>%
  mutate_at(c('braakstage1', 'braakstage2'), as.numeric) # convert chr to numeric


# average BRAAK stage columns
#* if between diagnosis range select most severe eg 1.5 = 2
Brain_pathology = Brain_pathology %>%
  mutate(BraakStage = braakstage1+braakstage2) %>%
  mutate(BraakStage = if_else(braakstage2 == 1 | braakstage2 == 2 | braakstage2 == 3 | braakstage2 == 4 | braakstage2 == 5 | braakstage2 == 6, BraakStage/2, braakstage1)) %>%
  mutate(BraakStage = ifelse(is.na(braakstage1) == FALSE & is.na(braakstage2) == TRUE, braakstage1, BraakStage)) %>%
  mutate(BraakStage = ceiling(BraakStage)) %>%
  select(-braakstage1, -braakstage2)


# Convert Thalphase into ordinal Predictor ---------------------------------

# Get Average of Thalphase value and convert to ordinal predictor
#* If between diagnosis ranges round up eg 1.5 = 2
Brain_pathology = Brain_pathology %>%
  separate(thalphase, into = c("thalphase1","thalphase2")) %>%
  mutate_at(c("thalphase1","thalphase2"), as.numeric) %>%
  mutate(thalphase = thalphase1+thalphase2) %>%
  mutate(thalphase = ifelse(is.na(thalphase2) == FALSE, thalphase/2, thalphase1)) %>%
  mutate(thalphase = ceiling(thalphase)) %>%
  select(-thalphase1, -thalphase2)


# Convert clinical diagnosis to a binary predictor for AD ---------------

# Identify clinical diagnosis
unique(Brain_pathology$Clinicaldiagnosis)

# convert to binary predictor
Brain_pathology = Brain_pathology %>%
  mutate(clinicalAD = NA) %>%
  mutate(ClinicalAD = ifelse(str_detect(Clinicaldiagnosis, "AD") |
                             str_detect(Clinicaldiagnosis, "dementia") |
                             str_detect(Clinicaldiagnosis, "Dementia") |
                             str_detect(Clinicaldiagnosis, "Alzheimer"), TRUE, FALSE))


# Tidy Brain weight -------------------------------------------------------

# Convert Brain weight to numeric predictor
Brain_pathology = Brain_pathology %>%
  mutate(brainweight = parse_number(brainweight, na = c("n/a", "na")))


# Convert CERAD score to ordinal predictor --------------------------------
#* If predictor was between two diagnosis, most severe diagnosis was chosen

Brain_pathology = Brain_pathology %>%
  mutate(CERAD_Ordinal = NA) %>%
  mutate(CERAD_Ordinal = ifelse(str_detect(cerad, "0") == TRUE, 0, CERAD_Ordinal)) %>%
  mutate(CERAD_Ordinal = ifelse(str_detect(cerad, "A") == TRUE, 1, CERAD_Ordinal)) %>%
  mutate(CERAD_Ordinal = ifelse(str_detect(cerad, "B") == TRUE, 2, CERAD_Ordinal)) %>%
  mutate(CERAD_Ordinal = ifelse(str_detect(cerad, "C") == TRUE, 3, CERAD_Ordinal)) %>%
  mutate(CERAD_Ordinal = ceiling(CERAD_Ordinal))


# Combine and convert both Neuropathical results into single categorical 
# predictor -------------------------------------------------------------------
#* Categories include 'AD', 'ADsymptoms', 'None' and NA

#* Identify unique phrases used in first neuropathical diagnosis
unique(Brain_pathology$Neuropathdiagnosis1)

# Convert into categories three categories
Brain_pathology = Brain_pathology %>%
  mutate(Neuropath1 = "") %>%
  mutate(Neuropath1 = ifelse(str_detect(Neuropathdiagnosis1, "incipient Alzheimer's") |
                                str_detect(Neuropathdiagnosis1, "mild Alzheimer") |
                                str_detect(Neuropathdiagnosis1, "ossible Alzheimer's") |
                                str_detect(Neuropathdiagnosis1, "Probable Alzheimer's") |
                                str_detect(Neuropathdiagnosis1, "mild AD") |
                                str_detect(Neuropathdiagnosis1, "probable Alzheimer's disease") |
                                str_detect(Neuropathdiagnosis1, "Incipient Alzheimer's disease") |
                                str_detect(Neuropathdiagnosis1, "early /incipient AD") |
                                str_detect(Neuropathdiagnosis1, "ARTAG") |
                                str_detect(Neuropathdiagnosis1, "DLB (Limbic; Braak 4)") |
                                str_detect(Neuropathdiagnosis1, "Possible Alzheimer's") |
                                str_detect(Neuropathdiagnosis1, "Alpha-synucleinopathy") |
                                str_detect(Neuropathdiagnosis1, "Alzheimer's disease (intermediate)") |
                                str_detect(Neuropathdiagnosis1, "early/incipient AD") |
                                str_detect(Neuropathdiagnosis1, "erebral amyloid") |
                                str_detect(Neuropathdiagnosis1, "ild tau")
                              == TRUE, "ADsymptoms", Neuropath1)) %>%
  mutate(Neuropath1 = ifelse(str_detect(Neuropathdiagnosis1, "Probable Alzheimer's disease") |
                                str_detect(Neuropathdiagnosis1, "moderate AD pathology") |
                                str_detect(Neuropathdiagnosis1, "Alzheimer's Disease") |
                                str_detect(Neuropathdiagnosis1, "robable AD") |
                                str_detect(Neuropathdiagnosis1, "Alzheimer's disease") |
                                str_detect(Neuropathdiagnosis1, "Posterior cortical atrophy variant of AD") |
                                str_detect(Neuropathdiagnosis1, "Moderate AD changes in temporal lobe")
                              == TRUE, "AD", Neuropath1)) %>%
  mutate(Neuropath1 = ifelse(str_detect(Neuropathdiagnosis1, "Ageing related changes") |
                                str_detect(Neuropathdiagnosis1, "Age related changes (mild)") |
                                str_detect(Neuropathdiagnosis1, "age changes only") |
                                str_detect(Neuropathdiagnosis1, "Normal for age")
                              == TRUE, "None", Neuropath1)) %>%
  mutate(Neuropath1 = ifelse(str_detect(Neuropathdiagnosis1, "STILL IN PROCESS") == TRUE, NA, Neuropath1))


#* Identify unique phrases used in second neuropathical diagnosis
unique(Brain_pathology$Neuropathdiagnosis2)

# Convert into categories three categories
Brain_pathology = Brain_pathology %>%
  mutate(Neuropath2 = "") %>%
  mutate(Neuropath2 = ifelse(str_detect(Neuropathdiagnosis2, "ossible Alzheimer") |
                                str_detect(Neuropathdiagnosis2, "ild Alzheimer") |
                                str_detect(Neuropathdiagnosis2, "mod CAA") |
                                str_detect(Neuropathdiagnosis2, "ild AD") |
                                str_detect(Neuropathdiagnosis2, "moderate CAA") |
                                str_detect(Neuropathdiagnosis2, "mild medial temporal lobe tau") |
                                str_detect(Neuropathdiagnosis2, "early/incipient Alzheimer's disease") |
                                str_detect(Neuropathdiagnosis2, "moderate tauopathy") |
                                str_detect(Neuropathdiagnosis2, "Hippocampal sclerosis") |
                                str_detect(Neuropathdiagnosis2, "v.mild AD") |
                                str_detect(Neuropathdiagnosis2, "AD changes") |
                                str_detect(Neuropathdiagnosis2, "ARTAG") |
                                str_detect(Neuropathdiagnosis2, "Low probability of AD") |
                                str_detect(Neuropathdiagnosis2, "CAA$")
                              == TRUE, "ADsymptoms", Neuropath1)) %>%
  mutate(Neuropath2 = ifelse(str_detect(Neuropathdiagnosis2, "AD. Mod") |
                                str_detect(Neuropathdiagnosis2, "moderate to severe SVD, very severe CAA") |
                                str_detect(Neuropathdiagnosis2, "moderate AD pathology") |
                                str_detect(Neuropathdiagnosis2, "^Alzheimer's disease") |
                                str_detect(Neuropathdiagnosis2, "evere CAA") |
                                str_detect(Neuropathdiagnosis2, "Intermediate AD pathology")
                              == TRUE, "AD", Neuropath1)) %>%
  mutate(Neuropath2 = ifelse(str_detect(Neuropathdiagnosis2, "age changes only") |
                                str_detect(Neuropathdiagnosis2, "Normal for age") |
                                str_detect(Neuropathdiagnosis2, "(age related)") |
                                str_detect(Neuropathdiagnosis2, "(Ageing)")
                             == TRUE, "None", Neuropath1))

# Get concensus from neuropathological diagnosis 1 and 2 
Brain_pathology = Brain_pathology %>%
  mutate(NeuropathDiagnosisConcensus = "") %>%
  mutate(NeuropathDiagnosisConcensus = if_else(Neuropath1 != Neuropath2, "conflict", NeuropathDiagnosisConcensus)) %>%
  mutate(NeuropathDiagnosisConcensus = if_else(Neuropath1 == Neuropath2, Neuropath1, NeuropathDiagnosisConcensus)) %>%
  mutate(NeuropathDiagnosisConcensus = if_else(is.na(Neuropath2) == TRUE, Neuropath1, NeuropathDiagnosisConcensus))

#remove columns
colnames(Brain_pathology)

Brain_pathology = Brain_pathology %>%
  select(-dateofdeath, -approxtimeofdeath, -dateofbrainremoval,
         -Neuropathdiagnosis1, -Neuropathdiagnosis2, -Clinicaldiagnosis,
         -cerad, -Neuropath1, -Neuropath2)


# convert Neuropathdiagnosis Concensus to binary predictor AD or NotAD
Brain_pathology = Brain_pathology  %>% 
  mutate(NeuropathDiagnosisConcensus = if_else(NeuropathDiagnosisConcensus == "ADsymptoms", "1", NeuropathDiagnosisConcensus)) %>%
  mutate(NeuropathDiagnosisConcensus = if_else(NeuropathDiagnosisConcensus == "conflict", "1", NeuropathDiagnosisConcensus)) %>%
  mutate(NeuropathDiagnosisConcensus = if_else(NeuropathDiagnosisConcensus == "", "0", NeuropathDiagnosisConcensus)) %>%
  mutate(NeuropathDiagnosisConcensus = if_else(is.na(NeuropathDiagnosisConcensus)==TRUE, "0", NeuropathDiagnosisConcensus)) %>%
  mutate(NeuropathDiagnosisConcensus = if_else(NeuropathDiagnosisConcensus == "None", '0', NeuropathDiagnosisConcensus)) %>%
  mutate(NeuropathDiagnosisConcensus = if_else(NeuropathDiagnosisConcensus== "AD", "1", NeuropathDiagnosisConcensus)) %>%
  mutate(NeuropathDiagnosisConcensus = na_if(NeuropathDiagnosisConcensus, "unknown")) %>%
  mutate(NeuropathDiagnosisConcensus = as_factor(NeuropathDiagnosisConcensus))


# Modify apoe_4 to detect homozygous e4 genotype --------------------------
Brain_pathology = Brain_pathology  %>% 
  mutate(apoe_e4 = ifelse(apoe == '4_4', 2, apoe_e4))


# Convert numeric factors to ordinal --------------------------------------
Brain_pathology = Brain_pathology %>%
  mutate(apoe_e4 = factor(apoe_e4, order = TRUE, levels = c(0,1,2))) %>%
  mutate(BraakStage = factor(BraakStage, order = TRUE, levels = c(0,1,2,3,4,5,6))) %>%
  mutate(thalphase = factor(thalphase, order = TRUE, levels = c(0,1,2,3,4,5))) %>%
  mutate(CERAD_Ordinal = factor(CERAD_Ordinal, order = TRUE, levels = c(0,1,2,3))) %>%
  mutate(NeuropathDiagnosisConcensus = factor(NeuropathDiagnosisConcensus, order = TRUE, levels = c(0,1))) %>%
  select(-clinicalAD)




# check whether co factors are represented (example below Apoe)
xtabs(~NeuropathDiagnosisConcensus + apoe, Brain_pathology)
xtabs(~NeuropathDiagnosisConcensus + apoe_e4, Brain_pathology)


# check whether how closely clinical diagnosis compare
xtabs(~NeuropathDiagnosisConcensus + ClinicalAD, Brain_pathology)
xtabs(~NeuropathDiagnosisConcensus + apoe_e4, Brain_pathology)


############################
# 4. Impute missing values #
############################
#* Imputing vissing values can increase the power of analysis. A single
#* imputation method will be used on predictors, where there is reasonable
#* reliability.

# Identify columns with missing data
na.patterns = naclus(select(Brain_pathology, colnames(Brain_pathology)[apply(Brain_pathology, 2, anyNA)]))

# Plot missing data fraction and association
naplot(na.patterns, 'na per var')
plot(na.patterns)

# impute missing values using single imputation method
w = transcan(~ ageatdeath + PMD + brainweight + apoe_e4 + BraakStage + thalphase +
               ClinicalAD + CERAD_Ordinal + NeuropathDiagnosisConcensus,
             imputed = TRUE, trantab=TRUE, data = Brain_pathology, pl = FALSE, pr = FALSE)


# See R^2 value to determine how reliable imputed values are for each predictor 
w$rsq

# Plot distribution
ggplot(w)

#Impute missing values
#* Based on fraction of missing values and imputation R2 value the following 
#* predictors will be imputed
attach(Brain_pathology)

Brain_pathology = Brain_pathology %>%
    mutate(apoe_e4 = impute(w, apoe_e4, data = Brain_pathology)) %>%
    mutate(BraakStage = impute(w, BraakStage, data = Brain_pathology)) %>%
    mutate(thalphase = impute(w, thalphase, data = Brain_pathology)) %>%
    mutate(CERAD_Ordinal = impute(w, CERAD_Ordinal, data = Brain_pathology)) %>%
    mutate(NeuropathDiagnosisConcensus = impute(w, NeuropathDiagnosisConcensus, data = Brain_pathology))


###################
# 5. Examine data #
###################

# Ranking predictors association to (Neuropathological/Clinical) AD diagnosis 
#* Somers Rank 
sDxy1 = with(Brain_pathology, rcorrcens(NeuropathDiagnosisConcensus ~  ageatdeath + PMD + brainweight + apoe_e4 + BraakStage + thalphase +
                   CERAD_Ordinal))

sDxy2 = with(Brain_pathology, rcorrcens(ClinicalAD ~  ageatdeath +PMD + brainweight + apoe_e4 + BraakStage + thalphase +
                                         CERAD_Ordinal))

# plot to figures together
par(mfrow=c(1,2))

# plot Somers Rank
plot(sDxy1, main = 'Somers Rank: Neuropathology Diagnosis', cex = 0.8 )
plot(sDxy2, main = 'Somers Rank: Clinical Diagnosis', cex = 0.8)


# return default plotting parameters
par(mfrow=c(1,1))



# redundancy analysis
redundancy = with(Brain_pathology, redun(~ ageatdeath + PMD + brainweight + apoe_e4 + BraakStage +
                     thalphase + ClinicalAD + CERAD_Ordinal +
                     NeuropathDiagnosisConcensus, r2 = 0.7, type = 'adjusted'))

redundancy$Out



##################################
# 6. Prepare Dataframe for PLINK #
##################################

# Convert all missing data to '-9
Brain_pathology_tidy = Brain_pathology %>%
  mutate(across(where(is.numeric), ~replace_na(.x, -9))) %>%
  mutate(across(where(is.character), ~replace_na(.x, '-9')))


Brain_pathology_tidy = Brain_pathology_tidy  %>%
  mutate(IID = ID,
         FID = IID) %>%
  select(FID, IID, case:thalphase, CERAD_Ordinal, everything(), -c(ID, apoe))

# without APOE e4 risk group
Brain_pathology_APOEStrat = Brain_pathology_tidy  %>%
  filter(apoe_e4 == "0")

###################
# 7. Export Files #
###################
#* export as textfile

# all data
write.table(Brain_pathology_tidy2,
            file = '../HLA_logitudinal_study/HLA_Plink_files/TidyPhenotype/BrainPathology_tidy.txt',
            quote = FALSE)


# Subset non APOE4 
write.table(Brain_pathology_APOEStrat,
            file = '../HLA_logitudinal_study/HLA_Plink_files/TidyPhenotype/Brain_pathology_APOEStrat.txt',
            quote = FALSE)

