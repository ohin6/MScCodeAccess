###############
# Description #
###############
#* Comparing HLA allele frequencies to reference population for quality control.
#* Reference populations are
      #* 1. "USA NMDP European Caucasian"
      #* 2. "USA NMDP Chinese"
      #* 3. "USA NMDP African American pop 2"
#* 
#* Alternative populations include: 
      #* 1. "USA NMDP Hispanic South or Central American"
      #* 2. "USA NMDP Japanese" 
      #* 3. "USA NMDP North American Amerindian"
      #* 4. "USA NMDP South Asian Indian"

#* Population used in study is from a Caucasian UK population and therefore it 
#* is expected to most closely follow HLA typing allele frequencies from a 
#* European Caucasian population.
#* 
#* This involves:
#* 1. Importing HLA allele population - formatted to MiDAS package requirements
#* 2. Compare HLA allele frequencies against reference population
#* 3. plot HLA frequency comparisons 
#* 
#* Requirements:
#* Need to have run MakeHLAMatrix.R script first to get corrected HLA typing of
#* population in a format usuable to MiDAS package
#* 
## Author: Dr. Owen Williams
##
## Date Created: 05-12-2022
##
## Email: owen.williams8@nhs.net

####################
# Install packages #
####################

library("tidyverse")
library("midasHLA")
library('stringr')

#########################
# Get working directory #
#########################
# get working directory
setwd(here::here())

#########################
# 1. Import HLA dataset #
#########################

# Import HLA genotypes
genotype = read.table('../HLA_logitudinal_study/HLA_Plink_files/updated_HLA_alleleTyping/MiDASGeno.txt')

# genotype = read.table('../Raw_data/MiDAS/MiDASGeno.txt')
colnames(genotype) = genotype[1,]
genotype = genotype[-1,]

##################################################################
# 2. Compare HLA allele frequencies against reference population #
##################################################################


# 2.1. Import reference panel HLA population frequencies -------------
freq_HLA <- getHlaFrequencies(hla_calls = genotype, compare = TRUE) %>%
  filter(Freq > 0.01)

# Identify reference populations
colnames(freq_HLA[4:10])


# 2.2 Wrangle data and select populations to compare freq against --------------
freq_HLA_long <- tidyr::gather(
  data = freq_HLA,
  key = "population",
  value = "freq",
  "Freq",
  "USA NMDP European Caucasian",
  "USA NMDP Chinese",
  "USA NMDP African American pop 2",
  factor_key = TRUE
)


# 2.3  Create function to compare against allele type---------------------------
#* Inputs i) HLA frequency dataset and reference pop, ii) HLA allele
#* output: plot
Compare_HLA_allele_freq = function(HLAfreq_data, HLA_allele){
  # filter allele
  HLAfreq_filter = HLAfreq_data %>%
    filter(grepl(paste0('^', HLA_allele), allele))
  
  # Create plot
  plot_HLAfreq =
  ggplot2::ggplot(data = HLAfreq_filter, ggplot2::aes(x = allele, y = freq, fill = population)) +
  ggplot2::geom_bar(
    stat = "identity",
    position = ggplot2::position_dodge(0.7),
    width = 0.7,
    colour = "black"
  ) +
  ggplot2::coord_flip() +
  ggplot2::scale_y_continuous(labels = formattable::percent) +
    ggplot2::theme(axis.text.y = element_text(size = 8, angle = 45,))
  return(plot_HLAfreq)
}

#########################
# 3. plot HLA frequency #
#########################

# Plot frequencies of specific alleles
plot1 = Compare_HLA_allele_freq(HLAfreq_data = freq_HLA_long, HLA_allele = 'A')
plot2 = Compare_HLA_allele_freq(HLAfreq_data = freq_HLA_long, HLA_allele = 'B')
plot3 = Compare_HLA_allele_freq(HLAfreq_data = freq_HLA_long, HLA_allele = 'C')
plot4 = Compare_HLA_allele_freq(HLAfreq_data = freq_HLA_long, HLA_allele = 'DRB1')

require(ggpubr)

x = ggarrange(plot1 + rremove("xylab"), plot2 + rremove("xylab"),
          plot3 + rremove("xylab"), plot4 + rremove("xylab"), 
          # + rremove("x.text"), 
          labels = c("a)", "b)", "c)", "d)"),
          ncol = 2, nrow = 2,
          align = 'v',
          common.legend = TRUE, legend = "top")


annotate_figure(x,
                bottom = text_grob("Frequency"),
                left = text_grob("HLA Alleles", rot = 90)
)


