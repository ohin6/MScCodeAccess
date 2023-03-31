###############
# Description #
###############
#* Script perform colocalisation between GWAS and eQTL data in the HLA region
#* using microglia single-cell sequencing data
#* 
#* GWAS data includes all studies which significant variant was identified which
#* include
#* 1. Bellenguez
#* 2. Kunkle
#* 3. Lambert
#* 4. Jansen
#* 5. Marioni
#* 6. Schwartzentruber
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
#* 4. Get GWAS study 
#* 5. Get LD matrix
#* 6. Compile datasets for coloc 
#* 7. Coloc for Single Causal variant
#* 8. Repeat for other GWAS studies

#* WARNING: These functions depend on LDlinkR API and can fail if server is
#* overloaded. If this happens repeat function until it works
#* 
#* WARNING: As this is working with large data some of these functions may take 
#* a long time populate.

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
library("coloc")
library("seqminer")
require(ggpubr)
require(LDlinkR)

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
# 3. Extract eQTL #
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

#* Examine HLA-DRB1 expression in microglia (repeat for others genes in subsequent analysis)
microglia = HLA_expression_tissue("microglia", DRB1)

###########################
# Convert to GRCh37 build #
###########################

# Convert column names to standardised format
sumStats = microglia[[2]] %>%
  rename(CHR = chromosome) %>%
  rename(SNP = variant) %>%
  rename(BP = position) %>%
  rename(P = pvalue)

# Convert to GRCh38 to GRCh37 using mungestats lift over
sumStats = MungeSumstats::liftover(sumstats_dt = sumStats,
                                   ref_genome = "hg38",
                                   convert_ref_genome = "hg19")


##################
# Get GWAS study #
##################

# Lambert_2013 ------------------------------------------------------------
# Get GWAS HLA region
Lambert_HLA_loci = data.table::fread('../Raw_data/Lambert_2013.tsv', header = T) %>%
  filter(str_detect(chromosome, '^6')) %>%
  filter(between(base_pair_location,32300000,32750000))

Lambert_HLA_loci = Lambert_HLA_loci %>%
  rename(CHR=chromosome,
         BP = base_pair_location,
         SNP= variant_id,
         BETA = beta,
         SE = standard_error,
         P = p_value)

Lambert_HLA_loci = Lambert_HLA_loci %>%
  mutate(P = as.numeric(P)) %>%
  filter(str_detect(CHR, '^6')) %>%
  filter(!str_detect(SNP, '^6'))
#################################################
# Filter both studies so they contain same SNPS #
#################################################
#* Include SNPs found in both studies based on genomic position. Also reduce
#* number of SNPs to < 1000 to get Linkage Disequillibrium

Lambert_HLA_loci2 <- Lambert_HLA_loci[order(Lambert_HLA_loci$BP),] %>%
  filter(BP %in% sumStats$BP) %>%
  dplyr::filter(P<0.00007)

sumStats =  sumStats[order(sumStats$BP),]  %>%
  filter(BP %in% Lambert_HLA_loci2$BP)

#################
# Get LD matrix #
#################

#LDlink token = b68c33126cbf
token = 'b68c33126cbf'

LD = LDlinkR::LDmatrix(snps = Lambert_HLA_loci2$SNP, pop = "EUR", r2d = "r2",
                       token = token, file = TRUE)


# filter SNP list from inputted GWAS file and remove any missing variants not
# found in imported LD matrix output
Lambert_HLA_loci2 = Lambert_HLA_loci2 %>%
  filter(SNP %in% as.vector(LD[1]$RS_number))


# repeat for eQTL data
sumStats =  sumStats[order(sumStats$BP),] %>%
  filter(BP %in% Lambert_HLA_loci2$BP)

# remove SNP names as column
LD.x = LD[-1]
# rownames = colnmaes
rownames(LD.x) = Lambert_HLA_loci2$BP
colnames(LD.x) = rownames(LD.x)
# Convert to matrix and remove NAs
LD.x = data.matrix(LD.x %>%
                     mutate(across(everything(), ~ replace_na(.x, 0))))



###############################
# Compile datasets for coloc #
##############################

# order SNPs accroding to position
sumStats = sumStats[order(sumStats$BP),]

# eQTL data set
D1 = list(beta = sumStats$beta,
          varbeta = sumStats$se*sumStats$se,
          SNP = sumStats$BP,
          position = sumStats$BP,
          type = 'quant',
          LD = LD.x,
          sdY = 1)

D2 = list(beta = Lambert_HLA_loci2$BETA,
          varbeta = Lambert_HLA_loci2$SE*Lambert_HLA_loci2$SE,
          SNP = Lambert_HLA_loci2$BP,
          position = Lambert_HLA_loci2$BP,
          type = 'cc',
          LD = LD.x,
          sdY = 1)

################################### 
# Coloc for Single Causal variant #
###################################
#* Perform colocalisation with the assumption there is a single causal variant
my.res <- coloc.abf(dataset1=D1,
                    dataset2=D2,
                    p1= 1e-4,
                    p2 = 1e-4,
                    p12 = 1e-5)

# Results
print(my.res)

subset(my.res$results,SNP.PP.H4>0.015)

o <- order(my.res$results$SNP.PP.H4,decreasing=TRUE)
cs <- cumsum(my.res$results$SNP.PP.H4[o])
w <- which(cs > 0.95)[1]
my.res$results[o,][1:w,]$snp

sensitivity(my.res,"H4 > 0.9")

#* Getting poor results as there is likely more than 1 causal variant

######################################
# Coloc for Multiple causal variants #
######################################
#* Perform colocalisation analysis using the SuSiE regression framework that allows
#* multiple causal variants 

# Run Susie on both datasets
zscore1 = sumStats$beta/sumStats$se
S1 = susieR::susie_rss(zscore1, data.matrix(LD.x), n =974)
summary(S1)


zscore2 = Lambert_HLA_loci2$BETA/Lambert_HLA_loci2$SE
S2 = susieR::susie_rss(zscore2, data.matrix(LD.x), n =974)
summary(S2)

# Perform colocalisation
x =

if(requireNamespace("susieR",quietly=TRUE)) {
  susie.res=coloc.susie(S2,S1)
  print(susie.res$summary)
}

write_csv(x, '/Users/owen/Desktop/microglia.csv')

################################################################################
############ Repeat for different GWAS data and gene profiles ##################
################################################################################
