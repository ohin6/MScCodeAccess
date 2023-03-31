###############
# Description #
###############
#* Fine-mapping causal variants using SuSiE finemapping package, from summary
#* statistics data, including importing linkage reference panel from LDlink
#* 
#* Steps
      #* 1. Step up LD link to extract LD matrix
      #* 2. Import and filter GWAS files
      #* 3. Create Function to automate LD download and SuSiE finemapping per study
      #* 4. Create Function to plot LD heatmap per study
      #* 5. Fine map GWAS files
      #* 6. Produce LDheat maps
      #* 7. Concatenate plots
  
#* WARNING: These functions depend on LDlinkR API and can fail if server is
#* overloaded. If this happens repeat function until it works
#* 
#* WARNING: As this is working with large data some of these functions may take 
#* a long time populate.

## Author: Dr. Owen Williams
##
## Date Created: 31-01-2023
##
## Email: owen.williams8@nhs.net

####################
# Install packages #
####################
require(data.table)
require(tidyverse)
require(susieR)
require(LDlinkR)
require(gaston)

#################
# 1. Set up LDlink #
#################

#LDlink token = b68c33126cbf
token = 'b68c33126cbf'

########################################
# 2. Import and filter GWAS summary stats #
########################################
#* Define region to less than 1000 SNPs so LD matrix can queried (LDmatrix limit
#* = 1000)

# Marioni_2018 ------------------------------------------------------------
# Get GWAS HLA region
Maroni_HLA_loci = fread('../Raw_data/Marioni_2019.txt', header = T) %>%
  filter(str_detect(CHR, '^6')) %>%
  filter(between(BP,32300000,32750000))

# Plot region to ensure area is captured correctly
plot(-log10(Maroni_HLA_loci$P)~Maroni_HLA_loci$BP)

# Reduce to < 1000 based on pvalue (to be able to get LD panel)
MaroniSNP_list = Maroni_HLA_loci %>%
  dplyr::filter(str_detect(Maroni_HLA_loci$SNP, 'rs')) %>%
  dplyr::filter(P<0.004)

# replot region containing SNPs to be fine-mapped
plot(-log10(MaroniSNP_list$P)~MaroniSNP_list$BP)


# Kunkle ------------------------------------------------------------------

# Get GWAS HLA region
Kunkle_loci = fread('../Raw_data/Kunkle.txt', header = T) %>%
  filter(str_detect(Chromosome, '^6')) %>%
  filter(between(Position,32300000,32750000)) %>%
  mutate(Pvalue = as.numeric(Pvalue))

# Standarise column names
colnames(Kunkle_loci) = c('CHR', 'BP', 'SNP', 'A1', 'A2', 'BETA', 'SE', 'P')

# Plot region to ensure area is captured correctly
plot(-log10(Kunkle_loci$P)~Kunkle_loci$BP)

# Reduce to < 1000 based on pvalue (to be able to get LD panel)
KunkleSNP_list = Kunkle_loci %>%
  dplyr::filter(str_detect(SNP, 'rs')) %>%
  dplyr::filter(P<0.00004)

# replot region containing SNPs to be fine-mapped
plot(-log10(MaroniSNP_list$P)~MaroniSNP_list$BP)


# Jansen ------------------------------------------------------------------

# Get GWAS HLA region
Jansen_loci = fread('../Raw_data/Jansenetal_2019.txt', header = T) %>%
  filter(str_detect(CHR, '^6')) %>%
  filter(between(BP,32300000,32750000)) %>%
  mutate(P = as.numeric(P))

# Plot region to ensure area is captured correctly
plot(-log10(Jansen_loci$P)~Jansen_loci$BP)

# Reduce to < 1000 based on pvalue (to be able to get LD panel)
JansenSNP_list = Jansen_loci %>%
  dplyr::filter(str_detect(Jansen_loci$SNP, 'rs')) %>%
  dplyr::filter(P<0.000005)

# replot region containing SNPs to be fine-mapped
plot(-log10(JansenSNP_list$P)~JansenSNP_list$BP)


# Schwartzentruber_2021 ---------------------------------------------------

# Get GWAS HLA region
Schwartzentruber_loci = fread('../Raw_data/Schwartzentruber_2021.tsv', header = T) %>%
  filter(str_detect(chromosome, '^6')) %>%
  filter(between(base_pair_location,32300000,32750000)) %>%
  mutate(p_value = as.numeric(p_value))

# standardise column names
colnames(Schwartzentruber_loci) = c('variant_id', 'P', 'CHR', 'BP', 'A1', 'A2',
                                    'MAF', 'BETA', 'SE', 'SNP', 'GWAS_BETA', 
                                    'GWAS_SE', 'GWAS_P', 'GWAX_UKBB_BETA',
                                    'GWAX_UKBB_SE', 'GWAX_UKBB_P', 'DIRECT', 'I2',
                                    'HET_P', 'INFO')

# Plot region to ensure area is captured correctly
plot(-log10(Schwartzentruber_loci$P)~Schwartzentruber_loci$BP)

# Reduce to < 1000 based on pvalue (to be able to get LD panel)
SchwartzentruberSNP_list = Schwartzentruber_loci %>%
  dplyr::filter(str_detect(Schwartzentruber_loci$SNP, 'rs')) %>%
  dplyr::filter(P<0.000000018)

# replot region containing SNPs to be fine-mapped
plot(-log10(SchwartzentruberSNP_list$P)~SchwartzentruberSNP_list$BP)



# Bellenguez_2022 ---------------------------------------------------------
#* Because in GRch38 build need to convert to GRCh37
Bellenguez_loci = sqldf::read.csv.sql(file = '../Raw_data/Bellenguez_2022_buildGRCh38.tsv',
                                      sql = 'select * from file where chromosome = 6',
                                      sep = '\t')

# data wrangling
Bellenguez_loci = Bellenguez_loci %>%
  dplyr::select(variant_id, chromosome, base_pair_location, p_value, effect_allele, other_allele, effect_allele_frequency, beta,
                standard_error) %>%
  mutate(chromosome = paste0('chr', chromosome))

# Standardize column names
colnames(Bellenguez_loci) = c('SNP', 'CHR', 'BP', 'P', 'A1', 'A2', 'MAF', 'BETA',
                              'SE')


# Convert build from hg38 to hg19
Bellenguez_loci <- MungeSumstats::liftover(sumstats_dt = Bellenguez_loci, 
                                           ref_genome = "hg38",
                                           convert_ref_genome = "hg19")

# Get HLA region
Bellenguez_loci = Bellenguez_loci %>%
  filter(between(BP,32300000,32750000))
  

# Plot region to ensure area is captured correctly
plot(-log10(Bellenguez_loci$P)~Bellenguez_loci$BP)

# Reduce to < 1000 based on pvalue (to be able to get LD panel)
BellenguezSNP_list = Bellenguez_loci %>%
  dplyr::filter(str_detect(Bellenguez_loci$SNP, 'rs')) %>%
  dplyr::filter(P<0.00000708)

# replot region containing SNPs to be fine-mapped
plot(-log10(BellenguezSNP_list$P)~BellenguezSNP_list$BP)



# Lambert_2013 ------------------------------------------------------------

Lambert_HLA_loci = fread('../Raw_data/Lambert_2013.tsv', header = T) %>%
  filter(str_detect(chromosome, '^6')) %>%
  filter(between(base_pair_location,32300000,32750000)) %>%
  select(-c(10:12)) %>%
  mutate(p_value = as.numeric(p_value))

# standardise column names
colnames(Lambert_HLA_loci) = c('CHR', 'BP', 'SNP', 'A1', 'A2', 'BETA', 'SE', 'P', 'MAF')

# Plot region to ensure area is captured correctly
plot(-log10(Lambert_HLA_loci$P)~Lambert_HLA_loci$BP)

# Reduce to < 1000 based on pvalue (to be able to get LD panel)
LambertSNP_list = Lambert_HLA_loci %>%
  dplyr::filter(str_detect(Lambert_HLA_loci$SNP, 'rs')) %>%
  dplyr::filter(P<0.000024)

# replot region containing SNPs to be fine-mapped
plot(-log10(LambertSNP_list$P)~LambertSNP_list$BP)


############################################
# 3. Develop function to automate Finemapping #
############################################
#* Create a function which automates SuSiE fine mapping using inputs of whole
#* region of interest (Locus which is around Chr6:32300000-32750000) and SNPs of
#* interest (SNP_list < 1000 SNPs).
#* 
#* Using SNPs of interest this function will download suitable reference LD panel
#* from LDLink using European population, calculate Z-score for each variant and 
#* perform SuSiE fine mapping to identify credible sets with 95% confidence of 
#* non zero effect.
#* 
#* Inputs:
        #* 1) Loci = Filtered GWAS file containing region of interest around HLA
        #*           region around Chr6:32300000-32750000
        #* 
        #* 2) SNP_list =  Filtered GWAS file containing region of interest
        #*                around HLA region which has less than < 1000 SNPs
        #*        

#* Outputs include:
        #* i) SuSiE_rss fit,   
        #* ii) highlighted SNPs plot,
        #* iii) table of credible sets
        #* iv) info of credible sets


get_cc = function(SNP_list, Locus){
#* Import LD matrix using LD link
LD = LDlinkR::LDmatrix(snps = SNP_list$SNP, pop = "EUR", r2d = "r2",
              token = token, file = TRUE)

# filter SNP list from inputted GWAS file and remove any missing variants 
SNP_list = SNP_list %>%
  filter(SNP %in% as.vector(LD[1]$RS_number))

# Convert LD into a matrix and convert NAs to 0 to make SuSiE compatabilble
LDmatrix = data.matrix(LD[-1] %>%
                         mutate(across(everything(), ~ replace_na(.x, 0))))

# specify how many samples included
n=nrow(LD)

# calculate z-score
zscore = SNP_list$BETA/SNP_list$SE

# Set random seed to allow for repeatability
set.seed(1)

# Perform SUSIE analysis using reference LD panel EUR
fit = susie_rss(zscore, sqrt(LDmatrix), n=n)

# Plot posterior inclusion probability
Plot_PIP = susie_plot(fit, y='PIP', add_bar = T,add_legend = T)

# identify location of causal variants
sets = fit$sets

# Plot
#* highlight based on credible sets
highlight = SNP_list %>%
  slice(c(unlist(fit$set$cs)))
# plot
plot_causal = ggplot(Locus, aes(x=BP,y=-log10(P))) +
  geom_point() +
  geom_label(
    data = highlight,
    nudge_x = 50000,
    aes(label = SNP)) +
  geom_point(data=highlight, 
             aes(BP,-log10(P)), 
             color='red',
             size=3) +
  theme_classic() +
  geom_hline(yintercept=5, linetype="dashed", 
             color = "red")

return(list(fit, plot_causal, highlight, sets))
}


#####################################
# 4. Develop function to LD heatmap #
#####################################
#* Using SNP_list of
LDheat = function(SNP_list, filename){
  
  # For computational effecencies Reduce dataframe size by choosing every third row
  SNP_list = SNP_list[seq(1, nrow(SNP_list), 3), ]
  
  # order by position descending order
  SNP_list = SNP_list[order(BP),]
  
  LD = LDlinkR::LDmatrix(snps = SNP_list$SNP, pop = "EUR", r2d = "r2",
                         token = token, file = FALSE)
  
  # filter SNP list from inputted GWAS file and remove any missing variants 
  SNP_list = SNP_list %>%
    filter(SNP %in% as.vector(LD[1]$RS_number))
  
  # Convert first row as rownames
  rownames(LD) = LD[,1]
  
  # Convert LD into a matrix and convert NAs to 0 
  LDmatrix = data.matrix(LD[-1] %>%
                           mutate(across(everything(), ~ replace_na(.x, 0))))
  
  return(LD.plot(LDmatrix, snp.positions = SNP_list$BP, pdf.file = paste0(filename,".pdf")))
}


#########################
# 5 Perform SuSiE model #
#########################
#* Perform finemapping on GWAS studies using 

# Marioni
SusieMarioni = get_cc(MaroniSNP_list, Maroni_HLA_loci)
SusieMarioni[[2]]
susie_plot(SusieMarioni[[1]], y='PIP', add_bar = T,add_legend = T)


# Kunkle
SusieKunkle = get_cc(KunkleSNP_list, Kunkle_loci)
SusieKunkle[[2]]

# Jansen
SusieJansen = get_cc(JansenSNP_list, Jansen_loci)
SusieJansen[[3]]

# Schwartzentruber
SusieSchwartzentruber = get_cc(SchwartzentruberSNP_list, Schwartzentruber_HLA_loci)
SusieSchwartzentruber[[2]]
susie_plot(SusieSchwartzentruber[[1]], y='PIP', add_bar = T,add_legend = T)

# Lambert_2013
SusieLambert = get_cc(LambertSNP_list, Lambert_HLA_loci)
SusieLambert[[2]]
susie_plot(SusieLambert[[1]], y='PIP', add_bar = T,add_legend = T)
SusieLambert[[4]]

# Bellenguez
SusieBellenguez = get_cc(BellenguezSNP_list, Bellenguez_loci)
SusieBellenguez[[2]]
susie_plot(SusieBellenguez[[1]], y='PIP', add_bar = T,add_legend = T)
SusieBellenguez[[4]]

###########################
# 6. Produce LD heat maps #
###########################

LDheat(MaroniSNP_list, 'LDMarioni_2018')
LDheat(KunkleSNP_list, 'LDKunkle')
LDheat(JansenSNP_list, 'Jansen')
LDheat(SchwartzentruberSNP_list, 'Schwartzentruber')
LDheat(BellenguezSNP_list, 'LDBellenguez')
LDheat(LambertSNP_list, 'LDlambert')


########################
# 7. Concatenate plots #
########################

# Manhatten plots -------------------------------------------------------

SusieBellenguez[[2]]
SusieJansen[[2]]
SusieKunkle[[2]]
SusieLambert[[2]]
SusieMarioni[[2]]
SusieSchwartzentruber[[2]]


# PIP ---------------------------------------------------------------------

susie_plot(SusieBellenguez[[1]], y='PIP', add_bar = T,add_legend = T)
susie_plot(SusieJansen[[1]], y='PIP', add_bar = T,add_legend = T)
susie_plot(SusieKunkle[[1]], y='PIP', add_bar = T,add_legend = T)
susie_plot(SusieLambert[[1]], y='PIP', add_bar = T,add_legend = T)
susie_plot(SusieMarioni[[1]], y='PIP', add_bar = T,add_legend = T)
susie_plot(SusieSchwartzentruber[[1]], y='PIP', add_bar = T,add_legend = T)
dev.off()



