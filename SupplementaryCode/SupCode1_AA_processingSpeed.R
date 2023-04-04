###############
# Description #
###############
#* 
#* Association analysis on cognitive decline predictor 'longitudingal Processing
#* Speed' (gsstd_lin) using data from the ACPRC cognitive decline cohort study.
#* Using a linear regression corrected for age and sex

#* Analysis was performed using the MiDAS package across different HLA 
#* groupings:
#* 1. Allele (low/high resolution)
#* 2. SuperType
#* 3. Functional epitope
#* 
#* Sample size = 1560
#*  
#* 
#* Steps:
#* 1. Import data
#* 2. Create custom function
#* 3. Association analysis
#* 
## Author: Dr. Owen Williams
##
## Date Created: 12-12-2022
##
## Email: owen.williams8@nhs.net


###########
# Summary #
###########

#* Strong association was observed between Supertype B62 and reduced processing 
#* speed, albeit just outside signiifcance threshold. See HLA_Supertypes table

####################
# Install packages #
####################

require(midasHLA)
require(tidyverse)
require(rms)

###################
# Import datasets #
###################

#  get HLA data
#* remove any white spaces

# Low HLA resolution (2 digits)
dat_HLA_LowRes <-
  readHlaCalls(
    file = '../HLA_logitudinal_study/HLA_Plink_files/updated_HLA_alleleTyping/MiDASGeno.txt',
    resolution = 2
  )

# High HLA resolution (4 digits)
dat_HLA_HighRes <-
  readHlaCalls(
    file = '../HLA_logitudinal_study/HLA_Plink_files/updated_HLA_alleleTyping/MiDASGeno.txt',
    resolution = 4
  )


# Phenotype data
cognitive = read.csv('/Users/owen/Library/Mobile Documents/com~apple~CloudDocs/MSc Project/R-scripts/Raw_data/CohortStudy/PhenotypicData/umlcha_sleepcog_04-2021x.csv') %>%
  select(FID, p1age, sex, gfstd_int:gvstd_lin) %>%
  mutate(across(everything(), ~replace(.x, .x == -9, NA))) %>%
  rename(ID = FID) %>%
  mutate(gfstd_int = as.numeric(gvstd_int))


###################
# Create function #
###################
#* Create function to run associations analysis with user input, including 
#* i) HLA genotype (matrix), ii) Phenotype dataframe, iii) MiDAS experiment, 
#* iv) Association model


# Create function
#* Note function uses i) dominant inheritance model, ii) 0.02-0.98 MAF cut off 
#* and iii) user selected p.adjust (recommend toggle between 'fdr' and 'bonferroni')

HLA_Associations = function(HLA_df, Pheno_df, experiment, Model, OmnibusTest, p.adjustTest){
  # Combine Geno and phenotype data
  MiDAS_Object = prepareMiDAS(
    hla_calls = HLA_df,
    colData = Pheno_df,
    experiment = experiment)
  
  # Create Model
  HLA_model = eval(parse(text=Model))
  
  # Get results
  HLA_results <- runMiDAS(
    object = HLA_model, 
    experiment = tail(experiment, n=1), # If more than one experiment choose last 
    inheritance_model = "dominant",
    omnibus = OmnibusTest,
    lower_frequency_cutoff = 0.02, 
    upper_frequency_cutoff = 0.98, 
    correction = p.adjustTest, 
    exponentiate = TRUE
  )
  return(HLA_results)
}


# Define HLA model
#* Do not change 'data' variable
model = 'glm(gsstd_lin ~ term + rcs(p1age,3) + sex,
                data = MiDAS_Object, family = gaussian())'



########################
# Association analysis #
########################


# Allele level at low resolution
HLA_LowRes = HLA_Associations(dat_HLA_LowRes, cognitive, "hla_alleles", model, FALSE, 'fdr')
HLA_LowRes

# Allele level at high resolution
HLA_HighRes = HLA_Associations(dat_HLA_HighRes, cognitive, "hla_alleles", model, FALSE, 'fdr')
HLA_HighRes

# Supertypes level
HLA_SuperTypes = HLA_Associations(dat_HLA_HighRes, cognitive, "hla_supertypes", model, FALSE, 'fdr')
HLA_SuperTypes
# B62 associated with AD

# Functional Epitope
NKcell = HLA_Associations(dat_HLA_HighRes, cognitive, c("hla_alleles", "hla_NK_ligands"), model, FALSE, 'fdr')
NKcell

write_csv(HLA_LowRes, '../../write up/Supplementary/SuppTable15_alleleLowRes.csv')
write_csv(HLA_HighRes, '../../write up/Supplementary/SuppTable15_alleleHighRes.csv')
write_csv(HLA_SuperTypes, '../../write up/Supplementary/SuppTable15_alleleSupertypes.csv')
write_csv(NKcell, '../../write up/Supplementary/SuppTable15_FunctionEpitope.csv')
