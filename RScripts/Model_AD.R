###############
# Description #
###############
#* The purpose of this script is developed an optimised logistic regression
#* prognostic model for predicting neuropathological AD diagnosis using data
#* ACPRC study.
#* 
#* Using previous AD classification performed in tidydata.R (n = 137) and
#* integrating cognitive decline APOE data information to inform model.
#* 
#* Steps:
#* 
#* 1. Import Data set
#* 2. Integrate data APOE and brain pathology results
#* 3. Tidy Data
#* 4. Subset data into Training and to be predicted data sets 
#* 5. Impute missing data - Model
#* 6. Identify predictors for Model
#* 7. Develop Model
#* 8. Internal Validation and Calibration
#* 9. Summarise Model
#* 10. Predict AD in external dataset
#* 
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
require(glmnet) # lasso regression
require(ROCit)


#####################
# 1. Import dataset #
#####################

# set current directory
setwd(here::here())

# Import dataset
Cognition <- read_csv("../Raw_data/CohortStudy/PhenotypicData/umlcha_sleepcog_04-2021x.csv") %>%
  dplyr::select(FID, sex, p1age, gfstd_int:gvstd_lin) %>%
  mutate(across(where(is.numeric), ~replace(.x, .x == -9, NA))) %>%
  mutate(across(where(is.character), ~replace(.x, .x== '-9', NA)))


# Remove rows where there are more than 4 missing data
cnt_na <- apply(Cognition, 1, function(z) sum(is.na(z)))
Cognition = Cognition[cnt_na < 4,]

##############################################
# 2. Integrate APOE and Brain phenotype data #
##############################################

# Import APOE data
APOE = read_excel("../Raw_data/CohortStudy/PhenotypicData/APOE_E4pres_geno.xlsx") %>%
  dplyr::select(panelid, apoe_e4, apoe)

# Import AD diagnosis
BrainPheno = read.table("../HLA_logitudinal_study/HLA_Plink_files/TidyPhenotype/BrainPathology_tidy.txt", header = T) %>%
  dplyr::select(FID, NeuropathDiagnosisConcensus)

# Join phenotype data
Cognition = Cognition %>%
  left_join(BrainPheno, c("FID" = "FID")) %>%
  left_join(APOE, c("FID" = "panelid"))

# Remove APOE and BrainPheno files to save memory
rm(APOE, BrainPheno)

################
# 3. Tidy data #
################

# Modify apoe_4 to detect homozygous e4 genotype --------------------------
Cognition = Cognition  %>%
  mutate(apoe_e4 = ifelse(apoe == '4_4', 2, apoe_e4)) %>%
  mutate(apoe_e4 = factor(apoe_e4, levels = c(0,1,2))) %>%
  dplyr::select(FID, apoe_e4, everything(), -apoe)

# modify colunm names
colnames(Cognition) = c('FID', 'APOE', 'Sex', 'Age', 'gfstd_int', 'gfstd_lin',
                        'gmstd_int', 'gmstd_lin', 'gsstd_int', 'gsstd_lin',
                        'gvstd_int', 'gvstd_lin', 'AD')


#########################################################
# 4. Subset data into modeling and predicting data sets #
#########################################################
#* Data which has AD diagnosis will be used for Modeling (model) while data
#* where AD diagnosis is absent will be Predicted (predict) 

model = subset(Cognition, !is.na(AD))
predict =subset(Cognition, is.na(AD))

###########################################
# 5. Impute missing values for model data #
###########################################
#* Imputing missing values can increase the power of analysis. A single
#* imputation method will be used on predictors, where there is reasonable
#* reliability.

attach(model)
# Create data distribution object
dd = datadist(model); options(datadist = 'dd')

# impute missing values using single imputation method
w = transcan(~ APOE + Sex + Age + gfstd_int + gfstd_lin + gmstd_int +
               gmstd_lin + gsstd_int + gsstd_lin + gvstd_int + gvstd_lin + AD,
             imputed = TRUE, trantab=TRUE, data = model, pl = FALSE, pr = FALSE)


# See R^2 value to determine how reliable imputed values are for each predictor 
w$rsq

#* Missing values for APOE will be imputed as the R2 value is high
model = model %>%
  mutate(APOE = impute(w, APOE, data = model))



###################################################
# 6. Identifying which predictors to use in model #
###################################################

# Convert APOE to binary predictor
#* Due to low sample size in 'model' dataset
model = model  %>%
  mutate(APOE = ifelse(APOE == 2, 1, APOE)) %>%
  mutate(APOE = factor(APOE, levels = c(0,1)))


#* Somers rank usedf to show which variables indivdually are most correlated to
#* response

sDxy = with(model,
             rcorrcens(AD ~ APOE + factor(Sex) + Age + gfstd_int + gfstd_lin + gmstd_int +
                         gmstd_lin + gsstd_int + gsstd_lin + gvstd_int + gvstd_lin))

plot(sDxy)

#* Redundancy analysis to determine redundant predictors (r2 at 0.4)
redundancy = with(model,
                  redun(~ factor(Sex) + Age + gfstd_int +
                          gfstd_lin + gmstd_int + gmstd_lin + gsstd_int +
                          gsstd_lin + gvstd_int + gvstd_lin,
                        r2 = 0.4, type = 'adjusted'))

redundancy$Out

# Lasso regression to identify most informative markers

# Predictors selected for lasso
lassoPred_comp = c('Age', 'gfstd_int', 'gfstd_lin', 'gmstd_int', 'gmstd_lin', 'gsstd_int',
                          'gsstd_lin', 'gvstd_int', 'gvstd_lin')

# Perform lasso regression
set.seed(1)
df = model %>% dplyr::select(dplyr::matches(lassoPred_comp))
lassoFit = cv.glmnet(as.matrix(df), as.matrix(AD), type.measure = "mse")
plot(lassoFit)

# Print coefficient values
lassoFit

# fit model using lasso regression calculated as lambda min or lambda to 1 SE
df2 = model %>% dplyr::select(matches(lassoPred_comp))
f.glmnet = glmnet(df2, AD, lambda = exp(-3))

# get coefficients of optimised model
Lassocoef = coef(f.glmnet, s = exp(-3))
Lassocoef = as.matrix(Lassocoef)
Lassocoef = tibble(predictor = row.names(Lassocoef), Lassocoef)
colnames(Lassocoef) = c("predictor", "coef")

# Show only included predictors
Lassocoef %>%
  filter(coef != 0)

#* Based off Lasso regression and Somers rank it is determined that i) gfstd_lin
#* ii) gmstd_lin, iii) gsstd_lin and iv) Sex are the most informative markers for
#* this model


#######################
# 7. Developing Model #
#######################
#* Develop model for predicting AD diagnosis 
f = lrm(factor(AD) ~ gfstd_lin + gmstd_lin + gsstd_lin + Sex, data=model)


# Assess Model
#* Anova
anova(f)

#########################################
# 8. Internal validation and calibration #
##########################################

#* Validate model
f = update (f, x=TRUE , y=TRUE)
v = validate (f, B=1000)

print(v, B=20, digits =3)

# Calcualte concordance index
Dxy = v[1,5]
c.index = Dxy*0.5+0.5

#* Calibrate Model
cal = calibrate(f, B=1000)
plot(cal)

######################
# 9. Summarise Model #
######################

# plot odds ratio
plot(summary(f), log =TRUE)

# Plot nomogram
plot(nomogram(f,
          fun=plogis , funlabel ="AD Risk Probability ",
          fun.at =c(.01 ,.05 ,.1 ,.25 ,.5 ,.75 ,.9 ,.95 ,.99 )))


####################
# 10. ROC analysis #
####################

model$prediction <- predict(f, type = "fitted")

# ROCit objects -----------------------------------------------------------

# Empirical
rocit_emp <- rocit(score = model$prediction,
                   class = model$AD, 
                   method = "emp")
# Binomial
rocit_bin <- rocit(score = model$prediction,
                   class = model$AD, 
                   method = "bin")

# Summary of ROC analysis -------------------------------------------------

summary(rocit_emp)
summary(rocit_bin)

# Confidence intervals 95%
set.seed(1)
ciROC_bin95 = ciROC(rocit_bin, 
                    level = 0.95, nboot = 1e4)
set.seed(1)
ciROC_emp95 = ciROC(rocit_emp, 
                    level = 0.95, nboot = 1e4)

# Area under the curve confidence interval
set.seed(1)
ciAUC_empboot = ciAUC(rocit_emp, 
                      nboot = 1e4)
print(ciAUC_empboot)

# Get cutoff values for ROC analysis
cutoff = data.frame(cbind(Cutoff=rocit_emp$Cutoff, 
                    TPR=rocit_emp$TPR, 
                    FPR=rocit_emp$FPR))

optimal = cutoff %>%
  filter(TPR > 0.63 & TPR < 0.7) %>%
  filter(FPR > 0.4 & FPR < 0.43)

print(optimal)

message('Optimal cutoff based on Youden index is: ', optimal[3,1])


# Plot ROC ----------------------------------------------------------------
#Set figure format
par(mfrow=c(1,2))

plot(rocit_emp, col = c(2,4), 
             legend = FALSE)
        lines(rocit_bin$TPR~rocit_bin$FPR, 
              col = 1, lwd = 2)
        lines(ciROC_emp95$LowerTPR~ciROC_emp95$FPR, 
              col = 2, lty = 2)
        lines(ciROC_emp95$UpperTPR~ciROC_emp95$FPR, 
              col = 2, lty = 2)
legend("bottomright", col = c(2,2,1,4),
       lty = c(1,2,1,2),
       lwd = c(2,1,2,2),
       c("Empirical ROC", "95% CI", "Binormal ROC", 'Chance line'))


# Compare with cumulative density function plot ---------------------------
#* cumulative density functions of postive and negative populations

kplot = ksplot(rocit_emp)

# summarise kplot
kplot$`KS stat`
kplot$`KS Cutoff` # optimal cut-off 

# Get optimal cut off
cutoff[which.min(abs(0.4846084-cutoff$Cutoff)),]
  

# postive to negative at optimal cutoff -----------------------------------

table(Actualvalue=model$AD,Predictedvalue=model$prediction>0.48)

dev.off()

######################################
# 11. Predict AD in external dataset #
######################################
#* Filter external dataset to include predictors used in model and drop rows
#* with missing data

predict = predict %>%
  dplyr::select(FID, gfstd_lin, gmstd_lin, gsstd_lin, Sex) %>%
  drop_na()

# Calculate predicted AD values and convert
#* This is given as probability where 1 AD and 0 is not AD
predict = predict %>%
  mutate(AD_diagnosis_prob = predict(f, data.frame(predict), type="fitted")) %>%
# Convert probabilities to binary response 1 = AD, 0 = not AD
  mutate(AD = ifelse(AD_diagnosis_prob >= 0.48, 1, 0))

# Export dataset
write.csv(predict, '../HLA_logitudinal_study/HLA_Plink_files/TidyPhenotype/Ext_predictedAD.csv')
write.csv(predict, '../../write up/Supplementary/Supp_Tables/SupTab1_Cohort.csv')

