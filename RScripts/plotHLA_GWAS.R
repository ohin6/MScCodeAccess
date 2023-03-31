###############
# Description #
###############
#* Creating Manhatten plots from the filtered summary statistics files created 
#* in GWAStoGenomicRangesChrm6.R.
#* 
#* This involves 
      #* 1. importing chrm6 genomic ranges summary stats files
      #* 2. Create Manhatten plots for all GWAS files
      #* 3. Create Manhatten plots around HLA region
#* 
#* previously filtered summary statistics files saved in directory:
#* '../Raw_data/Granges_GWASchrm6/' where there are 10 files
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
require(TxDb.Hsapiens.UCSC.hg19.knownGene)
require(karyoploteR)
require(sqldf)

##################################
# 1.Read in Gr_ranges if created #
##################################
#* Gr_ranges of chrm6 should have already been created using GWAStoGenomicRangesChrm6.R.
#* then these can be imported back in rather than needing to import large files
#* in above script.

# 1.1 Create function to convert datatable to genomic range ------------------

gRange = function(x){
  grange = makeGRangesFromDataFrame(x,
                                    seqnames.field="seqnames",
                                    start.field="start",
                                    end.field="end")
  grange$pval = x$pval
  names(grange) = x$rsid
  return(grange)
}

# 1.2 Run function on all files -----------------------------------------------

Backman_2021 = read_tsv('../Raw_data/Granges_GWASchrm6/Gr_Backman_2021_GRCh37.tsv')
Gr_Backman_2021 = gRange(Backman_2021)
rm(Backman_2021)

Bellenguez_2022 = read_tsv('../Raw_data/Granges_GWASchrm6/Gr_Bellenguez_2022_GRCh37.tsv')
Gr_Bellenguez_2022 = gRange(Bellenguez_2022)
rm(Bellenguez_2022)

Bi_2020 = read_tsv('../Raw_data/Granges_GWASchrm6/Gr_Bi_2020_GRCh37.tsv')
Gr_Bi_2020 = gRange(Bi_2020)
rm(Bi_2020)

Jansen_2019 = read_tsv('../Raw_data/Granges_GWASchrm6/Gr_Jansen_2019_GRCh37.tsv')
Gr_Jansen_2019 = gRange(Jansen_2019)
rm(Jansen_2019)

Kunkle_2019 = read_tsv('../Raw_data/Granges_GWASchrm6/Gr_Kunkle_2019_GRCh37.tsv')
Gr_Kunkle_2019 = gRange(Kunkle_2019)
rm(Kunkle_2019)

Lambert_2013 = read_tsv('../Raw_data/Granges_GWASchrm6/Gr_Lambert_2013_GRCh37.tsv')
Gr_Lambert_2013 = gRange(Lambert_2013)
rm(Lambert_2013)

Marioni_2018 = read_tsv('../Raw_data/Granges_GWASchrm6/Gr_Maroni_2019_GRCh37.tsv')
Gr_Marioni_2018 = gRange(Marioni_2018)
rm(Marioni_2018)

MorenoGrau_2019 = read_tsv('../Raw_data/Granges_GWASchrm6/Gr_MorenoGrau_2019_GRCh37.tsv')
Gr_MorenoGrau_2019 = gRange(MorenoGrau_2019)
rm(MorenoGrau_2019)

Schott_2016 = read_tsv('../Raw_data/Granges_GWASchrm6/Gr_Schott_2016_GRCh37.tsv')
Gr_Schott_2016 = gRange(Schott_2016)
rm(Schott_2016)

Schwartzentruber_2021 = read_tsv('../Raw_data/Granges_GWASchrm6/Gr_Schwartzentruber_2021_GRCh37.tsv')
Gr_Schwartzentruber_2021 = gRange(Schwartzentruber_2021)
rm(Schwartzentruber_2021)


########################################
# 2. Combine Manahatten plots at Chrm6 #
########################################
#* Create Manhatten plots using the GWAS information above. This process involves
#* i) Create colour range for GWAS data points, ii) creating reference genome 
#* hg19 build and iii) plotting data GWAS as Manhatten.
#* 
#* Due to the large figure size of combined Manhatten plots this will be split
#* into two figures kp1 and kp2

# 2.1. Add colour scheme to data points ----------------------------------------
#* Create function
points.col = function(grRange, startCol, endCol){
  pointsCol = -log10(grRange$pval)
  pointsCol = colByValue(pointsCol, colors=c(startCol, endCol))
  return(pointsCol)
}

#* GWAS data points with grey/red colour scheme
Lambert_2013_points.col = points.col(Gr_Lambert_2013, "#BBBBBB00", "red")
Marioni_2018_points.col = points.col(Gr_Marioni_2018, "#BBBBBB00", "red")
Kunkle_2019_points.col = points.col(Gr_Kunkle_2019, "#BBBBBB00", "red")
Bi_2020_points.col = points.col(Gr_Bi_2020, "#BBBBBB00", "red")
Backman_2021_points.col = points.col(Gr_Backman_2021, "#BBBBBB00", "red")

#* GWAS data points with grey/navy blue colour scheme
Schott_2016_points.col = points.col(Gr_Schott_2016, "#BBBBBB00", "navyblue")
Jansen_2019_points.col = points.col(Gr_Jansen_2019, "#BBBBBB00", "navyblue")
MorenoGrau_2019_points.col = points.col(Gr_MorenoGrau_2019, "#BBBBBB00", "navyblue")
Schwartzentruber_2021_points.col = points.col(Gr_Schwartzentruber_2021, "#BBBBBB00", "navyblue")
Bellenguez_2022_points.col = points.col(Gr_Bellenguez_2022, "#BBBBBB00", "navyblue")

# 2.2. Create reference plot 1-----------------------------------------------
kp <- plotKaryotype(plot.type=4, chromosomes = 'chr6', genome = 'hg19')
kpAddBaseNumbers(kp)

# 2.3. Add GWAS data plot 1-----------------------------------------------------------

kpAxis(kp, ymin=0, ymax= max(-log10(Gr_Lambert_2013$pval)), r0=autotrack(1,5), cex = 0.5)
kpAddLabels(kp, labels = "Lambert_2013", srt=90, pos=3, cex = .45, r0=autotrack(1,5))
kp = kpPlotManhattan(kp, Gr_Lambert_2013, points.col= Lambert_2013_points.col, r0=autotrack(1,5))

kpAxis(kp, ymin=0, ymax=max(-log10(Gr_Schott_2016$pval)), r0=autotrack(2,5), cex = 0.5)
kpAddLabels(kp, labels = "Schott_2016", srt=90, pos=3, cex = .45, r0=autotrack(2,5))
kp = kpPlotManhattan(kp, Gr_Schott_2016, points.col= Schott_2016_points.col, r0=autotrack(2,5))

kpAxis(kp, ymin=0, ymax=max(-log10(Gr_Marioni_2018$pval)), r0=autotrack(3,5), cex = 0.5)
kpAddLabels(kp, labels = "Marioni_2018", srt=90, pos=3, cex = .45, r0=autotrack(3,5))
kp = kpPlotManhattan(kp, Gr_Marioni_2018, points.col= Marioni_2018_points.col, r0=autotrack(3,5))

kpAxis(kp, ymin=0, ymax=max(-log10(Gr_Jansen_2019$pval)), r0=autotrack(4,5), cex = 0.5)
kpAddLabels(kp, labels = "Jansen_2019", srt=90, pos=3, cex = .45, r0=autotrack(4,5))
kp = kpPlotManhattan(kp, Gr_Jansen_2019, points.col= Jansen_2019_points.col, r0=autotrack(4,5))

kpAxis(kp, ymin=0, ymax=max(-log10(Gr_Kunkle_2019$pval)), r0=autotrack(5,5), cex = 0.5)
kpAddLabels(kp, labels = "Kunkle_2019", srt=90, pos=3, cex = .45, r0=autotrack(5,5))
kp = kpPlotManhattan(kp, Gr_Kunkle_2019, points.col = Kunkle_2019_points.col, r0=autotrack(5,5))


# 2.4. Add HLA region -------------------------------------------------------
#* Create Grange of region
hla_region = GRanges(seqnames = 'chr6',
                     ranges = IRanges(start = 29691116, end = 33054976))
#* Add to track
kpRect(kp, data=hla_region, y0=0, y1=1, col=NA, border="red", lwd=3)


# 2.5. Create reference plot 2 (repeat steps 2-4) -----------------------------
kp2 <- plotKaryotype(plot.type=4, chromosomes = 'chr6', genome = 'hg19')
kpAddBaseNumbers(kp2)

#* Add GWAS data plot 2 --------------------------------------------------------

kpAxis(kp2, ymin=0, ymax=max(-log10(Gr_MorenoGrau_2019$pval)), r0=autotrack(1,5), cex = 0.5)
kpAddLabels(kp2, labels = "MorenoGrau_2019", srt=90, pos=3, cex = .45, r0=autotrack(1,5))
kp2 = kpPlotManhattan(kp2, Gr_MorenoGrau_2019, points.col = MorenoGrau_2019_points.col, r0 = autotrack(1,5))

kpAxis(kp2, ymin=0, ymax=max(-log10(Gr_Bi_2020$pval)), r0=autotrack(2,5), cex = 0.5)
kpAddLabels(kp2, labels = "Bi_2020", srt=90, pos=3, cex = .45, r0=autotrack(2,5))
kp2 = kpPlotManhattan(kp2, Gr_Bi_2020, points.col=Bi_2020_points.col, r0 = autotrack(2,5))

kpAxis(kp2, ymin=0, ymax=max(-log10(Gr_Schwartzentruber_2021$pval)), r0=autotrack(3,5), cex = 0.5)
kpAddLabels(kp2, labels = "Schwartzentruber_2021", srt=90, pos=3, cex = .45, r0 = autotrack(3,5))
kp2 = kpPlotManhattan(kp2, Gr_Schwartzentruber_2021, points.col = Schwartzentruber_2021_points.col, r0=autotrack(3,5))

kpAxis(kp2, ymin=0, ymax=max(-log10(Gr_Backman_2021$pval)), r0=autotrack(4,5), cex = 0.5)
kpAddLabels(kp2, labels = "Backman_2021", srt=90, pos=3, cex = .45, r0=autotrack(4,5))
kp2 = kpPlotManhattan(kp2, Gr_Backman_2021, points.col = Backman_2021_points.col, r0=autotrack(4,5))

kpAxis(kp2, ymin=0, ymax=max(-log10(Gr_Bellenguez_2022$pval)), r0=autotrack(5,5), cex = 0.5)
kpAddLabels(kp2, labels = "Bellenguez_2022", srt=90, pos=3, cex = .45, r0=autotrack(5,5))
kp2 = kpPlotManhattan(kp2, Gr_Bellenguez_2022, points.col = Bellenguez_2022_points.col, r0=autotrack(5,5))

#* Add HLA region
kpRect(kp, data=hla_region, y0=0, y1=1, col=NA, border="red", lwd=3)


######################################
# 3. Combine plots around HLA region #
######################################

# 3.1. Set new colour scheme for data points -----------------------------------

#* GWAS data points with grey/red colour scheme
Lambert_2013_points.col = points.col(Gr_Lambert_2013, "#BBBBBB00", "red")
Marioni_2018_points.col = points.col(Gr_Marioni_2018, "#BBBBBB00", "navyblue")
Jansen_2019_points.col = points.col(Gr_Jansen_2019, "#BBBBBB00", "red")
Kunkle_2019_points.col = points.col(Gr_Kunkle_2019, "#BBBBBB00", "navyblue")
Schwartzentruber_2021_points.col = points.col(Gr_Schwartzentruber_2021, "#BBBBBB00", "red")
Bellenguez_2022_points.col = points.col(Gr_Bellenguez_2022, "#BBBBBB00", "navyblue")

# 3.2. Get SNP Peaks ----------------------------------------------------

#* Create function for identifying peaks
topSNP = function(Gr_range){
  # identify top SNP between two regions
  SNP1 = Gr_range[Gr_range@ranges@start>32000000 & Gr_range@ranges@start<33500000]
  SNP2 = Gr_range[Gr_range@ranges@start>40000000 & Gr_range@ranges@start<42000000]
  topSNPs = c(SNP1[which.min(SNP1$pval)],SNP2[which.min(SNP2$pval)])
  # Get y value (-log10(pval))
  topSNPs$y = -log10(topSNPs$pval)
  # remove top SNPs outside -log10 p-value threshold
  if(!is.null(which(topSNPs$y< -log10(5e-8))) == TRUE){
    topSNPs[which(topSNPs$y< -log10(5e-8))] = NULL
  }
  # Give Scale to y where max = 10
  maxY = max(topSNPs$y)
  topSNPs$y = topSNPs$y/maxY * 10
  return(topSNPs)
}

# Identify peaks on Granges

topSNP_Lambert_2013 = topSNP(Gr_Lambert_2013)
topSNP_Marioni_2018 = topSNP(Gr_Marioni_2018)
topSNP_Jansen_2019 = topSNP(Gr_Jansen_2019)
topSNP_Kunkle_2019 = topSNP(Gr_Kunkle_2019)
topSNP_Schwartzentruber_2021 = topSNP(Gr_Schwartzentruber_2021)
topSNP_Bellenguez_2022 = topSNP(Gr_Bellenguez_2022)


# 3.3 Set Region -----------------------------------------------------------
kp = plotKaryotype(plot.type=4, zoom="chr6:28000000-43000000")
kpAddBaseNumbers(kp, add.units = TRUE, cex=0.5, tick.dist = 1e6)

# 3.4. Highlight HLA region -------------------------------------------------
# Create Grange of region
hla_region = GRanges(seqnames = 'chr6',
                     ranges = IRanges(start = 29691116, end = 33054976))
# Add to track
kpRect(kp, data=hla_region, y0=0, y1=1, col=NA, border="red", lwd=3)

# 3.5. Create curated gene track --------------------------------------------
genes.data <- makeGenesDataFromTxDb(txdb = TxDb.Hsapiens.UCSC.hg19.knownGene, karyoplot = kp)
genes.data <- addGeneNames(genes.data)
genes.data <- mergeTranscripts(genes.data)
kpPlotGenes(kp, data=genes.data, add.transcript.names = FALSE, r1=0.25, cex=0.4, gene.name.position = "left", gene.name.cex = 0.4)
kpAddLabels(kp, labels = "Gene Track", srt=90, pos=3, r0=autotrack(1,4), cex = 0.7)

# 3.6. Add GWAS data plot 1 ----------------------------------------------------
# Add GWAS data
kpAxis(kp, ymin=0, ymax=max(-log10(topSNP_Jansen_2019$pval)), r0=autotrack(4,4), cex = 0.7)
kpAddLabels(kp, labels = "Jansen_2019", srt=90, cex = .5, pos=3, r0=autotrack(4,4))
kp = kpPlotManhattan(kp, Gr_Jansen_2019, points.col=Jansen_2019_points.col, r0=autotrack(4,4), highlight = topSNP_Jansen_2019)
kpText(kp, data = topSNP_Jansen_2019, labels = names(topSNP_Jansen_2019), ymax=10, pos=4, cex= 0.75, col="black", r0=autotrack(4,4))
kpPoints(kp, data = topSNP_Jansen_2019, pch=1, cex=1, col="black", lwd=2, ymax=10, r0=autotrack(4,4))

kpAxis(kp, ymin=0, ymax=max(-log10(topSNP_Marioni_2018$pval)), r0=autotrack(3,4), cex = 0.7)
kpAddLabels(kp, labels = "Marioni_2018", srt=90, cex = .5, pos=3, r0=autotrack(3,4))
kp = kpPlotManhattan(kp, Gr_Marioni_2018, points.col=Marioni_2018_points.col, r0=autotrack(3,4), highlight = topSNP_Marioni_2018)
kpText(kp, data = topSNP_Marioni_2018, labels = names(topSNP_Marioni_2018), ymax=10, pos=4, cex= 0.75, col="black", r0=autotrack(3,4))
kpPoints(kp, data = topSNP_Marioni_2018, pch=1, cex=1, col="black", lwd=2, ymax=10, r0=autotrack(3,4))

kp = kpPlotManhattan(kp, Gr_Lambert_2013, points.col=Lambert_2013_points.col, r0=autotrack(2,4), highlight = topSNP_Lambert_2013)
kpAxis(kp, ymin=0, ymax=max(-log10(topSNP_Lambert_2013$pval)), r0=autotrack(2,4), cex = 0.7)
kpAddLabels(kp, labels = "Lambert_2013", srt=90, pos=3, cex = .5, r0=autotrack(2,4))
kpText(kp, data = topSNP_Lambert_2013, labels = names(topSNP_Lambert_2013), ymax=10, pos=4, cex= 0.75, col="black", r0=autotrack(2,4))
kpPoints(kp, data = topSNP_Lambert_2013, pch=1, cex=1, col="black", lwd=2, ymax=10, r0=autotrack(2,4))

# 3.7 plot 2 ---------------------------------------------------------------
#* Specify Region
kp = plotKaryotype(plot.type=4, zoom="chr6:28000000-43000000")
kpAddBaseNumbers(kp, add.units = TRUE, cex=0.5, tick.dist = 1e6)
#* Add gene track
kpPlotGenes(kp, data=genes.data, add.transcript.names = FALSE, r1=0.25, cex=0.4, gene.name.position = "left", gene.name.cex = 0.4)
kpAddLabels(kp, labels = "Gene Track", srt=90, pos=3, r0=autotrack(1,4), cex = 0.7)
#* Add to HLA region
kpRect(kp, data=hla_region, y0=0, y1=1, col=NA, border="red", lwd=3)
#* Add GWAS data
kpAxis(kp, ymin=0, ymax=max(-log10(topSNP_Bellenguez_2022$pval)), r0=autotrack(4,4), cex = 0.7)
kpAddLabels(kp, labels = "Bellenguez_2022", srt=90, cex = .5, pos=3, r0=autotrack(4,4))
names(Gr_Bellenguez_2022) = NULL #change row names due to missing rsids
kp = kpPlotManhattan(kp, Gr_Bellenguez_2022, points.col = Bellenguez_2022_points.col, r0=autotrack(4,4), highlight = topSNP_Bellenguez_2022)
kpText(kp, data = topSNP_Bellenguez_2022, labels = names(topSNP_Bellenguez_2022), ymax=10, pos=4, cex= 0.75, col="black", r0=autotrack(4,4))
kpPoints(kp, data = topSNP_Bellenguez_2022, pch=1, cex=1, col="black", lwd=2, ymax=10, r0=autotrack(4,4))

kpAxis(kp, ymin=0, ymax=max(-log10(topSNP_Schwartzentruber_2021$pval)), r0=autotrack(3,4), cex = 0.7)
kpAddLabels(kp, labels = "Schwartzentruber_2021", srt=90, cex = .5, pos=3, r0=autotrack(3,4))
names(Gr_Schwartzentruber_2021) = NULL #change row names due to missing rsids
kp = kpPlotManhattan(kp, Gr_Schwartzentruber_2021, points.col = Schwartzentruber_2021_points.col, r0=autotrack(3,4), highlight = topSNP_Schwartzentruber_2021)
kpText(kp, data = topSNP_Schwartzentruber_2021, labels = names(topSNP_Schwartzentruber_2021), ymax=10, pos=4, cex= 0.75, col="black", r0=autotrack(3,4))
kpPoints(kp, data = topSNP_Schwartzentruber_2021, pch=1, cex=1, col="black", lwd=2, ymax=10, r0=autotrack(3,4))

kpAxis(kp, ymin=0, ymax=max(-log10(topSNP_Kunkle_2019$pval)), r0=autotrack(2,4), cex = 0.7)
kpAddLabels(kp, labels = "Kunkle_2019", srt=90, pos=3, cex = .5, r0=autotrack(2,4))
names(Gr_Kunkle_2019) = NULL #change row names due to missing rsids
kp = kpPlotManhattan(kp, Gr_Kunkle_2019, points.col=Kunkle_2019_points.col, r0=autotrack(2,4), highlight = topSNP_Kunkle_2019)
kpText(kp, data = topSNP_Kunkle_2019, labels = names(topSNP_Kunkle_2019), ymax=10, pos=4, cex= 0.75, col="black", r0=autotrack(2,4))
kpPoints(kp, data = topSNP_Kunkle_2019, pch=1, cex=1, col="black", lwd=2, ymax=10, r0=autotrack(2,4))



