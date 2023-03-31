###############
# Description #
###############
#* HLA association analysis of the cohort study using the MiDAS package.
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
dat_pheno = read.table('../HLA_logitudinal_study/HLA_Plink_files/TidyPhenotype/BrainPathology_tidy.txt', sep = '\t', header = T) %>%
  select(FID, ageatdeath:NeuropathDiagnosisConcensus) %>%
  rename(ID = FID) %>%
  mutate(across(everything(), ~replace(.x, .x=='-9', NA)))


###################
# Create function #
###################
#* Create function to run associations analysis with user input, including 
#* i) HLA genotype (matrix), ii) Phenotype dataframe, iii) MiDAS experiment, 
#* iv) Association model

# Define HLA model
#* Do not change 'data' variable
model = 'glm(NeuropathDiagnosisConcensus ~ term + rcs(ageatdeath,3),
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
HLA_LowRes = HLA_Associations(dat_HLA_LowRes, dat_pheno, "hla_alleles", model, FALSE, 'fdr')
HLA_LowRes

# Locus level at high resolution
HLA_HighRes = HLA_Associations(dat_HLA_HighRes, dat_pheno, "hla_alleles", model, FALSE, 'fdr')
HLA_HighRes

# Supertypes level
HLA_SuperTypes = HLA_Associations(dat_HLA_HighRes, dat_pheno, "hla_supertypes", model, FALSE, 'fdr')
HLA_SuperTypes
# B62 associated with AD



# amino acid level
HLA_aa = HLA_Associations(dat_HLA_HighRes, dat_pheno, "hla_aa", model, TRUE, 'fdr')
HLA_aa

#* Investigate effect estimates distributed for most significant ammino acid B_9:
MiDAS_Object_AA = prepareMiDAS(
  hla_calls = dat_HLA_HighRes,
  colData = dat_pheno,
  experiment = "hla_aa")

model_AA = glm(NeuropathDiagnosisConcensus ~ term + rcs(ageatdeath,3),
                data = MiDAS_Object_AA, family = binomial())

HLA_AA_B_9_results <- runMiDAS(
  model_AA,
  experiment = "hla_aa",
  inheritance_model = "dominant",
  omnibus_groups_filter = "B_9",
  lower_frequency_cutoff = 0.02,
  upper_frequency_cutoff = 0.98,
  correction = "fdr",
  exponentiate = TRUE
)

HLA_AA_B_9_results

#* Investigate effect estimates distributed for most significant ammino acid B_156:
HLA_AA_B_156_results <- runMiDAS(
  model_AA,
  experiment = "hla_aa",
  inheritance_model = "dominant",
  omnibus_groups_filter = "B_156",
  lower_frequency_cutoff = 0.02,
  upper_frequency_cutoff = 0.98,
  correction = "fdr",
  exponentiate = TRUE
)

HLA_AA_B_156_results

# Is there evidence to suggest AD is associated with NK cell interactions

NKcell = HLA_Associations(dat_HLA_HighRes, dat_pheno, c("hla_alleles", "hla_NK_ligands"), model, FALSE, 'fdr')
NKcell


#* Although not sigificant in this dataset HLA-Bw4 (on HLA-B) carriers are 
#* associated with increased disease risk of AD. Further KIR typing would be 
#* required to confirm