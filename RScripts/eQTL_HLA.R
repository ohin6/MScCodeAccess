###############
# Description #
###############
#* Script using Tabix to extract single cell sequencing eQTL data from the eQTL
#* catalogue online database.
#* 
#* The expression of HLA-DRB1 which has been linked to AD via the neuroinflammation
#* pathways has been was assessed across numerous tissue types across defined
#* HLA locus.
#* 
#* This achieved by creating a function which links to the eQTL catalogue and allows 
#* user to select study and tissue type to extract relevant information from.

#* Steps
#* 1. Set up link with online eQTL catalog
#* 2. Define QTL region
#* 3. Extract eQTL for different tissue types
#* 4. Create Manhatten plots for each tissue
#* 5. Get top SNP

## Author: Dr. Owen Williams
##
## Date Created: 21-11-2022
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

########################
# 2. Define QTL region #
########################
#* Region of interest is the HLA QTL on chromosome 6. The GTEx and Youngs_2019 
#* uses the GRCh38 build so this will need to be converted to GRCh37 build
#* afterwards. The region of interest in GRCh37 is Chr6:32300000-32750000 which 
#* when converted to GRCh38 is Chr6:32332223-32782223

# Define region of interest in GRCh38 build
region = "6:32332223-32782223"

#identify tissues in different studies
GTEx_studies = dplyr::filter(tabix_paths, study == "GTEx")
Young_studies = dplyr::filter(tabix_paths, study == "Young_2019")
Schmiedel_studies = dplyr::filter(tabix_paths, study == "Schmiedel_2018")

##############################################
# 3. Extract eQTL for different tissue types #
##############################################

# 3.1 Create function -----------------------------------------------------

#* Create function plot expression results in APOE region in user selected
#* tissue type (see above to see available tissues)
HLA_expression_tissue = function(tissue, STUDY){
  studies_df = dplyr::filter(tabix_paths, study == STUDY, tissue_label == tissue)
  #Extract column names from first file
  column_names = colnames(readr::read_tsv(studies_df$ftp_path, n_max = 1))
  #Import summary statistics - selected gene = HLA-DRB1
  summary_stats1 = import_eQTLCatalogue(studies_df$ftp_path, region, selected_gene_id = "ENSG00000196126", column_names) 
  
  # Convert column names to standardised format
  summary_stats1 = summary_stats1 %>%
    rename(CHR = chromosome) %>%
    rename(SNP = variant) %>%
    rename(BP = position) %>%
    rename(P = pvalue)
  
  # Convert to GRCh38 to GRCh37 using mungestats lift over
  summary_stats1 = MungeSumstats::liftover(sumstats_dt = summary_stats1,
                                     ref_genome = "hg38",
                                     convert_ref_genome = "hg19")
  # plot
  plot1 = ggplot(summary_stats1, aes(x = BP, y = -log10(P))) +
    labs(title = tissue) +
    geom_point() +
    theme_bw() +
    ylim(0,max(-log(summary_stats1$P,10)))
  list = list(plot1, summary_stats1)
  return(list)
}

#############################################
# 4. Create Manhatten plots for each tissue #
#############################################

# 3.2 Run function on specific tissues ---------------------------------------
# to be used in figure
Microglia = HLA_expression_tissue("microglia", 'Young_2019')
CD4 = HLA_expression_tissue("CD4+ T cell", 'Schmiedel_2018')
hippocampus = HLA_expression_tissue("brain (hippocampus)", 'GTEx')
cerebellum = HLA_expression_tissue("brain (cerebellum)", 'GTEx')
adiposeVisceral = HLA_expression_tissue("adipose (visceral)", 'GTEx')
hypothalamus = HLA_expression_tissue("brain (hypothalamus)", 'GTEx')
liver = HLA_expression_tissue("liver", 'GTEx')
frontal_cortex = HLA_expression_tissue("brain (DLPFC)", 'GTEx')

# Additional gene expression tissues
brain_cingulate = HLA_expression_tissue("brain (anterior cingulate cortex)", 'GTEx')
brain_putamen = HLA_expression_tissue("brain (putamen)", 'GTEx')
brain_substantiaNigra = HLA_expression_tissue("brain (substantia nigra)", 'GTEx')
brain_amygdala = HLA_expression_tissue("brain (amygdala)", 'GTEx')
brain_caudate = HLA_expression_tissue("brain (caudate)", 'GTEx')
brain_cortex = HLA_expression_tissue("brain (cortex)", 'GTEx')
brain_nucleusAccumbens = HLA_expression_tissue("brain (nucleus accumbens)", 'GTEx')
brain_spinalCord = HLA_expression_tissue("brain (spinal cord)", 'GTEx')
LCL = HLA_expression_tissue("LCL", 'GTEx')
Blood = HLA_expression_tissue("blood", 'GTEx')
adipose = HLA_expression_tissue("adipose", 'GTEx')
Kidney = HLA_expression_tissue("kidney (cortex)", 'GTEx')


# View plots --------------------------------------------------------------

Microglia[[1]] 
CD4[[1]]
hippocampus[[1]]
frontal_cortex[[1]]
brain_cingulate[[1]]
cerebellum[[1]]
hypothalamus[[1]]
brain_putamen[[1]]
brain_substantiaNigra[[1]]
brain_amygdala[[1]]
brain_caudate[[1]]
brain_cortex[[1]]
brain_hippocampus[[1]]
brain_nucleusAccumbens[[1]]
brain_spinalCord[[1]]
LCL[[1]]
Blood[[1]]
adiposeVisceral[[1]]
adipose[[1]]
Kidney[[1]]
liver[[1]]


# Create plot -------------------------------------------------------------

ggarrange(hippocampus[[1]], cerebellum[[1]], hypothalamus[[1]],
          frontal_cortex[[1]], liver[[1]], adiposeVisceral[[1]],
          Microglia[[1]], CD4[[1]],
          ncol = 4, nrow = 2)

#######################
# 5. Get leading SNPs #
#######################

leadSNPs = hippocampus[[2]] %>%
  dplyr::filter(P == min(P)) %>%
  mutate(Tissue = 'hippocampus')
leadSNPs = rbind(leadSNPs,
                 cerebellum[[2]] %>%
                   dplyr::filter(P == min(P))%>%
                   mutate(Tissue = 'cerebellum'),
                 hypothalamus[[2]] %>%
                   dplyr::filter(P == min(P))%>%
                   mutate(Tissue = 'hypothalamus'),
                 frontal_cortex[[2]] %>%
                   dplyr::filter(P == min(P))%>%
                   mutate(Tissue = 'frontal_cortex'),
                 liver[[2]] %>%
                   dplyr::filter(P == min(P))%>%
                   mutate(Tissue = 'liver'),
                 adiposeVisceral[[2]] %>%
                   dplyr::filter(P == min(P)) %>%
                   mutate(Tissue = 'adiposeVisceral'),
                 Microglia[[2]] %>%
                   dplyr::filter(P == min(P))%>%
                   mutate(Tissue = 'Microglia'),
                 CD4[[2]] %>%
                   dplyr::filter(P == min(P))%>%
                   mutate(Tissue = 'CD4+ T cell'))
                   
write_csv(leadSNPs, '/Users/owen/Desktop/leadSNPs.csv')
