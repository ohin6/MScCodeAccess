###############
# Description #
###############
#* From the downloaded GWAS summary statistics file this script i) converts to 
#* hg.19 Genome build, ii) Filters chromosome 6, and iii) converts region into a
#* GenomicRange file type.
#* 
#* The purpose of this is to store relevant information from the large GWAS
#* summary statistics files into smaller more manageable and accessible files
#* for later analysis
#* 
#* GWAS summary statistic files were previously downloaded from gwas catalog: 
#* (https://www.ebi.ac.uk/gwas/efotraits/MONDO_0004975) and saved to directory
#* '../Raw_data/'. In total there were 10 GWAS summary stats files
#* 
#* 
## Author: Dr. Owen Williams
##
## Date Created: 21-11-2022
##
## Email: owen.williams8@nhs.net

####################
# Install packages #
####################
require(tidyverse)
require(MungeSumstats)
require(GenomicRanges)

##########################################################
#  Download GWAS Data sets and convert to Genomic Ranges #
##########################################################
# Summary Statistics previously downloaded 

# 1.  Marioni_2018 ---------------------------------------------------------
# Read in file with SQL query
#* import chrm6 only
#* add sep = /t to covert to tsv
Marioni_2018 = sqldf::read.csv.sql(file = '../Raw_data/Marioni_2019.txt',
                                   sql = 'select * from file where CHR = 6',
                                   sep = ' ')

# wrangle data
Marioni_2018 = Marioni_2018 %>%
  dplyr::select(SNP, CHR, BP, P) %>%
  mutate(CHR = paste0('chr', CHR))

# Change col names
colnames(Marioni_2018) = c('SNP','seqnames', 'ranges', 'pval')

# Create GenomicRanges file
Gr_Marioni_2018 = makeGRangesFromDataFrame(Marioni_2018,
                                           seqnames.field=c("CHR"),
                                           start.field="BP",
                                           end.field=c("BP"))

# Add p-values
Gr_Marioni_2018$pval = Marioni_2018$P
# rsIDs
names(Gr_Marioni_2018) = Marioni_2018$SNP
# Remove dataframe to save memory
rm(Marioni_2018)

# write Grange to file
write.table(as.data.frame(Gr_Marioni_2018),
            file="../Raw_data/Granges_GWASchrm6/Gr_Maroni_2018_GRCh37.tsv", 
            quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t", na="")

# 2. Schwartzentruber_2021 ------------------------------------------------

# Import Chrm 6
Schwartzentruber_2021 = sqldf::read.csv.sql(file = '../Raw_data/Schwartzentruber_2021.tsv',
                                            sql = 'select * from file where chromosome = 6',
                                            sep = '\t')

# Data wrangling
Schwartzentruber_2021 = Schwartzentruber_2021 %>%
  mutate(chromosome = paste0('chr',chromosome)) %>%
  dplyr::select(variant_id, chromosome, base_pair_location, p_value)

# Standardize column names
colnames(Schwartzentruber_2021) = c('SNP', 'CHR', 'BP', 'P')

# Create GenomicRanges file
Gr_Schwartzentruber_2021 = makeGRangesFromDataFrame(Schwartzentruber_2021,
                                                    seqnames.field=c("CHR"),
                                                    start.field="BP",
                                                    end.field="BP")

# Add pvalue
Gr_Schwartzentruber_2021$pval = Schwartzentruber_2021$P
# rsIDs
names(Gr_Schwartzentruber_2021) = Gr_Schwartzentruber_2021$rsid

# Remove dataframe to save memory
rm(Schwartzentruber_2021)

# write Grange to file
write.table(as.data.frame(Gr_Schwartzentruber_2021),
            file="../Raw_data/Granges_GWASchrm6/Gr_Schwartzentruber_2021_GRCh37.tsv", 
            quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t", na="")

# 3. Kunkle_2019 ----------------------------------------------------------

# Import Chrm 6
Kunkle_2019 = sqldf::read.csv.sql(file = '../Raw_data/Kunkle.txt',
                                  sql = 'select * from file where chromosome = 6',
                                  sep = ' ')

# Data wrangling
Kunkle_2019 = Kunkle_2019 %>%
  dplyr::select(MarkerName, Chromosome, Position, Pvalue) %>%
  mutate(Chromosome = paste0('chr', Chromosome))

# Standardize column names
colnames(Kunkle_2019) = c('SNP', 'CHR', 'BP', 'P')

# Create GenomicRanges file
Gr_Kunkle_2019 = makeGRangesFromDataFrame(Kunkle_2019,
                                          ignore.strand=TRUE,
                                          seqnames.field=c("CHR"),
                                          start.field="BP",
                                          end.field="BP",
                                          starts.in.df.are.0based = FALSE)
# add p-values
Gr_Kunkle_2019$pval = Kunkle_2019$P
# rsIDs
names(Gr_Kunkle_2019) = Kunkle_2019$SNP


# Remove dataframe to save memory
rm(Kunkle_2019)

# write Grange to file
write.table(as.data.frame(Gr_Kunkle_2019),
            file="../Raw_data/Granges_GWASchrm6/Gr_Kunkle_2019_GRCh37.tsv", 
            quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t", na="")

# 4. Moreno-Grau ----------------------------------------------------------

# Import Chrm 6
MorenoGrau_2019 = sqldf::read.csv.sql(file = '../Raw_data/Moreno_Grau.txt',
                                      sql = 'select * from file where CHR = 6',
                                      sep = '\t')

# data wrangling
MorenoGrau_2019 =MorenoGrau_2019 %>%
  dplyr::select(rsID, CHR, Position, P) %>%
  mutate(CHR = paste0('chr', CHR))


# Standardize column names
colnames(MorenoGrau_2019) = c('SNP', 'CHR', 'BP', 'P')

# Create GenomicRanges file
Gr_MorenoGrau_2019 = makeGRangesFromDataFrame(MorenoGrau_2019,
                                              seqnames.field=c("CHR"),
                                              start.field="BP",
                                              end.field="BP")

# add p-values
Gr_MorenoGrau_2019$pval = MorenoGrau_2019$P
# add rsids
names(Gr_MorenoGrau_2019) = MorenoGrau_2019$SNP

# Remove dataframe to save memory
rm(MorenoGrau_2019)

# write Grange to file
write.table(as.data.frame(Gr_MorenoGrau_2019),
            file="../Raw_data/Granges_GWASchrm6/Gr_MorenoGrau_2019_GRCh37.tsv", 
            quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t", na="")

# 5. Jansen ---------------------------------------------------------------

Jansen_2019 = sqldf::read.csv.sql(file = '../Raw_data/Jansenetal_2019.txt',
                                  sql = 'select * from file where CHR = 6',
                                  sep = '\t')

# data wrangling
Jansen_2019 =Jansen_2019 %>%
  dplyr::select(SNP, CHR, BP, P) %>%
  mutate(CHR = paste0('chr', CHR))

# Standardize column names
colnames(Jansen_2019) = c('SNP', 'CHR', 'BP', 'P')

# Create GenomicRanges file
Gr_Jansen_2019 = makeGRangesFromDataFrame(Jansen_2019,
                                          seqnames.field=c("CHR"),
                                          start.field="BP",
                                          end.field="BP")

# add p-values
Gr_Jansen_2019$pval = Jansen_2019$P
# add rsids
names(Gr_Jansen_2019)= Jansen_2019$SNP

# Remove dataframe to save memory
rm(Jansen_2019)

# write Grange to file
write.table(as.data.frame(Gr_Jansen_2019),
            file="../Raw_data/Granges_GWASchrm6/Gr_Jansen_2019_GRCh37.tsv", 
            quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t", na="")

# 6. Bellenguez_2022_buildGRCh38 ------------------------------------------

Bellenguez_2022 = sqldf::read.csv.sql(file = '../Raw_data/Bellenguez_2022_buildGRCh38.tsv',
                                      sql = 'select * from file where chromosome = 6',
                                      sep = '\t')

# data wrangling
Bellenguez_2022 = Bellenguez_2022 %>%
  dplyr::select(variant_id, chromosome, base_pair_location, p_value) %>%
  mutate(chromosome = paste0('chr', chromosome))

# Standardize column names
colnames(Bellenguez_2022) = c('SNP', 'CHR', 'BP', 'P')


# Convert build from hg38 to hg19
Bellenguez_2022 <- MungeSumstats::liftover(sumstats_dt = Bellenguez_2022, 
                                           ref_genome = "hg38",
                                           convert_ref_genome = "hg19")


# Create GenomicRanges file
Gr_Bellenguez_2022 = makeGRangesFromDataFrame(Bellenguez_2022,
                                              seqnames.field=c("chromosome"),
                                              start.field="BP",
                                              end.field="BP")

# add p-values
Gr_Bellenguez_2022$pval = Bellenguez_2022$P
# add rsid
names(Gr_Bellenguez_2022) = Bellenguez_2022$SNP

# Remove dataframe to save memory
rm(Bellenguez_2022)

# write Grange to file
write.table(as.data.frame(Gr_Bellenguez_2022),
            file="../Raw_data/Granges_GWASchrm6/Gr_Bellenguez_2022_GRCh37.tsv", 
            quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t", na="")


# 7. Backman_2021_buildGRCh38 --------------------------------------------

# Import data at chrm6
Backman_2021 = sqldf::read.csv.sql(file = '../Raw_data/Backman_2021_GRCh38.tsv',
                                   sql = 'select * from file where chromosome = 6',
                                   sep = '\t')
# Data wrangling
Backman_2021 = Backman_2021 %>%
  dplyr::select(Name, chromosome, base_pair_location, p_value) %>%
  mutate(chromosome = paste0('chr', chromosome))

# Standardise column names
colnames(Backman_2021) = c('SNP', 'CHR', 'BP', 'P')

# Convert build from hg38 to hg19
Backman_2021 <- MungeSumstats::liftover(sumstats_dt = Backman_2021, 
                                        ref_genome = "hg38",
                                        convert_ref_genome = "hg19")

# change chromosome column to work with karyploteR package
Backman_2021 = Backman_2021 %>%
  mutate(CHR = paste0('chr', CHR))

# Create GenomicRanges file
Gr_Backman_2021 = makeGRangesFromDataFrame(Backman_2021,
                                           seqnames.field=c("CHR"),
                                           start.field="BP",
                                           end.field="BP")

# add p-values
Gr_Backman_2021$pval = Backman_2021$P
# Add rsids
names(Gr_Backman_2021)= Backman_2021$SNP

# Remove dataframe to save memory
rm(Backman_2021)

# write Grange to file
write.table(as.data.frame(Gr_Backman_2021),
            file="../Raw_data/Granges_GWASchrm6/Gr_Backman_2021_GRCh37.tsv", 
            quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t", na="")

# 8. BI_2020 --------------------------------------------------------------

# Read in file with SQL query
#* import chrm6 only
#* add sep = /t to covert to tsv
Bi_2020 = sqldf::read.csv.sql(file = '../Raw_data/Bi_2020.tsv',
                              sql = 'select * from file where chromosome = 6',
                              sep = '\t')

# wrangle data
Bi_2020 = Bi_2020 %>%
  select(variant_id, chromosome, base_pair_location, p_value, p_value_SPACox_NoSPA) %>%
  mutate(chromosome = paste0('chr', chromosome))

colnames(Bi_2020) = c('SNP', 'CHR', 'BP', 'P', 'P2')

# Create GenomicRanges file
Gr_Bi_2020 = makeGRangesFromDataFrame(Bi_2020,
                                      seqnames.field=c("CHR"),
                                      start.field="BP",
                                      end.field="BP")

# add p-values
Gr_Bi_2020$pval = Bi_2020$P
# Add rsid
names(Gr_Bi_2020) = Bi_2020$SNP

# Remove dataframe to save memory
rm(Bi_2020)

# write Grange to file
write.table(as.data.frame(Gr_Bi_2020),
            file="../Raw_data/Granges_GWASchrm6/Gr_Bi_2020_GRCh37.tsv", 
            quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t", na="")

# 9. Schott_2016 ----------------------------------------------------------

# Read in file with SQL query
#* import chrm6 only
#* add sep = ' ' to covert to tsv
Schott_2016 = sqldf::read.csv.sql(file = '../Raw_data/Schott_2016.tsv',
                                  sql = 'select * from file where CHR = 6',
                                  sep = ' ')

# wrangle data
Schott_2016 = Schott_2016 %>%
  dplyr::select(rsid, CHR, position, frequentist_add_pvalue) %>%
  mutate(CHR = paste0('chr', CHR))

# Standardise column names
colnames(Schott_2016) = c('SNP', 'CHR', 'BP', 'P')

# Create GenomicRanges file
Gr_Schott_2016 = makeGRangesFromDataFrame(Schott_2016,
                                          seqnames.field=c("CHR"),
                                          start.field="BP",
                                          end.field="BP")

# add p-values
Gr_Schott_2016$pval = Schott_2016$P
# add p-values
names(Gr_Schott_201) = Schott_2016$SNP


# Remove dataframe to save memory
rm(Schott_2016)

# write Grange to file
write.table(as.data.frame(Gr_Schott_2016),
            file="../Raw_data/Granges_GWASchrm6/Gr_Schott_2016_GRCh37.tsv", 
            quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t", na="")


# 10.Lambert_2013 ---------------------------------------------------------

# Read in file with SQL query
#* import chrm6 only
#* add sep = /t to covert to tsv
Lambert_2013 = sqldf::read.csv.sql(file = '../Raw_data/Lambert_2013.tsv',
                                   sql = 'select * from file where chromosome = 6',
                                   sep = '\t')

# wrangle data
Lambert_2013 = Lambert_2013 %>%
  select(variant_id, chromosome, base_pair_location, p_value) %>%
  mutate(chromosome = paste0('chr', chromosome))

# standardise column names
colnames(Lambert_2013) = c('SNP', 'CHR', 'BP', 'P')

# Create GenomicRanges file
Gr_Lambert_2013 = makeGRangesFromDataFrame(Lambert_2013,
                                           seqnames.field=c("CHR"),
                                           start.field="BP",
                                           end.field="BP")

# add p-values
Gr_Lambert_2013$pval = Lambert_2013$P
# add rsids
names(Gr_Lambert_2013) = Lambert_2013$SNP


# Remove dataframe to save memory
rm(Lambert_2013)

# write Grange to file
write.table(as.data.frame(Gr_Lambert_2013), file="../Raw_data/Granges_GWASchrm6/Gr_Lambert_2013_GRCh37.tsv", 
            quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t", na="")
