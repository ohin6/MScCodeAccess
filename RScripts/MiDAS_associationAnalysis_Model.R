###############
# Description #
###############
#* HLA association analysis using Model assisted AD diagnosis. Associations 
#* performed using MiDAS package
#* 
#* This script will examine associations at 
#* 
#* 1.  2 x Locus level
#* low resolution (HLA 2 digit)
#* High resolution (HLA 4 digit)
#* 2. Supertypes
#* 3. NK HLA interactions 
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

####################
# Install packages #
####################

require(midasHLA)
require(tidyverse)
require(rms)
require(here)

######################
# 1. Import datasets #
######################
# Set working directory
setwd(here::here())

# Import predicted phenotype data
Pred_pheno = read_csv('../HLA_logitudinal_study/HLA_Plink_files/TidyPhenotype/Ext_predictedAD.csv') %>%
  mutate(AD = factor(AD, levels= c(0,1))) %>%
  dplyr::rename(ID = FID) %>% 
  dplyr::select(ID:AD)

summary(Pred_pheno$AD)

# Import HLA data
#* Low resolution (2 digits)
dat_HLA_LowRes <-
  readHlaCalls(
    file = '../HLA_logitudinal_study/HLA_Plink_files/updated_HLA_alleleTyping/MiDASGeno.txt',
    resolution = 2
  )

# High resolution (4 digits)
dat_HLA_HighRes <-
  readHlaCalls(
    file = '../HLA_logitudinal_study/HLA_Plink_files/updated_HLA_alleleTyping/MiDASGeno.txt',
    resolution = 4
  )


###################
# Create function #
###################
#* Create function to run associations analysis with user input, including 
#* i) HLA genotype (matrix), ii) Phenotype dataframe, iii) MiDAS experiment, 
#* iv) Association model

# Define HLA model
#* Do not change 'data' variable
model = 'glm(AD ~ term,
                data = MiDAS_Object, family = binomial())'


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


########################
# Association analysis #
########################

# Locus level at low resolution
HLA_LowRes = HLA_Associations(dat_HLA_LowRes, Pred_pheno, "hla_alleles", model, FALSE, 'fdr')
HLA_LowRes
write_csv(HLA_LowRes, '../../write up/Supplementary/Pred_HLA_LowRes.csv')

# Locus level at high resolution
HLA_HighRes = HLA_Associations(dat_HLA_HighRes, Pred_pheno, "hla_alleles", model, FALSE, 'fdr')
HLA_HighRes
write_csv(HLA_HighRes, '../../write up/Supplementary/Pred_HLA_HighRes.csv')

# Supertypes level
HLA_SuperTypes = HLA_Associations(dat_HLA_HighRes, Pred_pheno, "hla_supertypes", model, FALSE, 'fdr')
HLA_SuperTypes
write_csv(HLA_SuperTypes, '../../write up/Supplementary/Pred_HLA_SuperTypes.csv')

# amino acid level
HLA_aa = HLA_Associations(dat_HLA_HighRes, Pred_pheno, "hla_aa", model, TRUE, 'fdr')
HLA_aa
write_csv(HLA_aa, '../../write up/Supplementary/Pred_HLA_aa.csv')

# Is there evidence to suggest AD is associated with NK cell interactions
NKcell = HLA_Associations(dat_HLA_HighRes, Pred_pheno, c("hla_alleles", "hla_NK_ligands"), model, FALSE, 'fdr')
NKcell
write_csv(NKcell, '../../write up/Supplementary/Pred_NKcell.csv')

