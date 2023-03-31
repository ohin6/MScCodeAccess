###############
# Description #
###############
#* Overlay fine-mapped variants against gene annotation data extracted from the
#* Notts et al., 2019 paper. Gene annotation data downloaded from the UCSC 
#* website: https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr6%3A32300000%2D32750000&hgsid=1574969533_7NnXzfZPW44XY95Et2DsG7vhcqor
#* 
#* This includes the (1) promoter, (2) enhancer and (3) superEnhancer of four
#* tissue types:

        #* 1. PU.1 (microglia)
        #* 2. LHx2 (hypothalic networks)
        #* 3. NeuN (Neuronal Nucleur antigen)
        #* 4. Oligodendrocyte (neuroglia)

#* as well as levels of Microglia
        #* 1. ATAC 
        #* 2. H3K4me3
        #* 3. chromatin arrangements 

#* Steps
#* 1. Import genome annotation files
#* 2. Create tracks 
#* 3. Import Reference Genes
#* 3. Create chromatin interaction track
#* 4. Create Plots

## Author: Dr. Owen Williams
##
## Date Created: 21-02-2023
##
## Email: owen.williams8@nhs.net

############
# Packages #
############

require(Gviz)
require(GenomicRanges)
require(tidyverse)
require(GenomicInteractions)

setwd(here::here())


##########################################
# 1. Import genome annotation from files #
##########################################
#* Previously downloaded to file. See description to find source URL

# Import datasets ---------------------------------------------------------
# Import LHX2
LHx2En = read.csv(file = 'LHx2enhncers.csv') %>% distinct() 
LHx2prom = read.csv(file = 'LHx2Promoters.csv') %>% distinct() 
LHX2SuperEn = read.csv(file = 'LHx2Superenhncers.csv') %>% distinct() 

# Import PU1
PU1En = read.csv(file = 'PU1Enhancers.csv') %>% distinct() 
PU1prom = read.csv(file = 'PU1Promoters.csv') %>% distinct() 
PU1SuperEn = read.csv(file = 'PU1Superenhncers.csv') %>% distinct() 

# Import Neuon
neuonEn = read.csv(file = 'neuonEnhancers.csv') %>% distinct() 
neuonProm = read.csv(file = 'neuonPromoters.csv') %>% distinct() 
neuonSuperEn = read.csv(file = 'neuonSuperenhncers.csv') %>% distinct() 

# Oligodendrocyte
olig2En = read.csv(file = 'olig2Enhancers.csv') %>% distinct() 
olig2prom = read.csv(file = 'olig2Promoters.csv') %>% distinct() 
olig2SuperEn = read.csv(file = 'olig2Superenhncers.csv') %>% distinct() 

# microglia superenhancer
MicroSup = read.csv(file = 'microgliaSuperEnhancer.csv') %>% distinct() 
# Microglia ATAC Pooled
MicroATAC = read_csv(file = 'microgliaAtac.csv') %>% distinct() %>%
  filter(!str_detect(chr, '#bedGraph'))
# PU1 H3K4 
PU1H3K4 = read_csv(file = 'PU1H3K4Me3.csv') %>% distinct() %>%
  filter(!str_detect(chr, '#bedGraph'))

####################
# 2. Create Tracks #
####################
#* Tracks for Enhancer (en), Promoter (prom) and Super Enhancer (SupEn)

# Tracks LHX2 -------------------------------------------------------------


LHx2En_Trk = AnnotationTrack(GRanges(seqnames = LHx2En$X.chrom,
                        ranges = IRanges(start = LHx2En$chromStart,
                                         end = LHx2En$chromEnd)),
                name = "LHX2 Enhancer")

LHx2prom_Trk = AnnotationTrack(GRanges(seqnames = LHx2prom$X.chrom,
                                     ranges = IRanges(start = LHx2prom$chromStart,
                                                      end = LHx2prom$chromEnd)),
                             name = "LHx2promoter")

LHx2SupEn_Trk = AnnotationTrack(GRanges(seqnames = LHX2SuperEn$X.chrom,
                                       ranges = IRanges(start = LHX2SuperEn$chromStart,
                                                        end = LHX2SuperEn$chromEnd)),
                               name = "LHX2 SuperEnhancer")


# Track PU1 ---------------------------------------------------------------


PU1En_Trk = AnnotationTrack(GRanges(seqnames = PU1En$X.chrom,
                                     ranges = IRanges(start = PU1En$chromStart,
                                                      end = PU1En$chromEnd)),
                             name = "PU.1 Enhancer")

PU1prom_Trk = AnnotationTrack(GRanges(seqnames = PU1prom$X.chrom,
                                       ranges = IRanges(start = PU1prom$chromStart,
                                                        end = PU1prom$chromEnd)),
                               name = "PU.1 promoter")

PU1SupEn_Trk = AnnotationTrack(GRanges(seqnames = PU1SuperEn$X.chrom,
                                        ranges = IRanges(start = PU1SuperEn$chromStart,
                                                         end = PU1SuperEn$chromEnd)),
                                name = "PU.1 SuperEnhancer")


# Track Oligodendrocyte ---------------------------------------------------


olig2En_Trk = AnnotationTrack(GRanges(seqnames = olig2En$X.chrom,
                                    ranges = IRanges(start = olig2En$chromStart,
                                                     end = olig2En$chromEnd)),
                            name = "Olig Enhancer")

olig2prom_Trk = AnnotationTrack(GRanges(seqnames = olig2prom$X.chrom,
                                      ranges = IRanges(start = olig2prom$chromStart,
                                                       end = olig2prom$chromEnd)),
                              name = "Olig promoter")

olig2SupEn_Trk = AnnotationTrack(GRanges(seqnames = olig2SuperEn$X.chrom,
                                       ranges = IRanges(start = olig2SuperEn$chromStart,
                                                        end = olig2SuperEn$chromEnd)),
                               name = "Olig SuperEnhancer")


# # Track neuron ----------------------------------------------------------


neuonEn_Trk = AnnotationTrack(GRanges(seqnames = neuonEn$X.chrom,
                                    ranges = IRanges(start = neuonEn$chromStart,
                                                     end = neuonEn$chromEnd)),
                            name = "NeuN Enhancer")

neuonprom_Trk = AnnotationTrack(GRanges(seqnames = neuonProm$X.chrom,
                                      ranges = IRanges(start = neuonProm$chromStart,
                                                       end = neuonProm$chromEnd)),
                              name = "NeuN promoter")

neuonSupEn_Trk = AnnotationTrack(GRanges(seqnames = neuonSuperEn$X.chrom,
                                       ranges = IRanges(start = neuonSuperEn$chromStart,
                                                        end = neuonSuperEn$chromEnd)),
                               name = "NeuN SuperEnhancer")



# Track Microglia -------------------------------------------------------


MicroSup_Trk = AnnotationTrack(GRanges(seqnames = MicroSup$X.chrom,
                                       ranges = IRanges(start = MicroSup$chromStart,
                                                        end = MicroSup$chromEnd)),
                               name = "Microglia Super Enhancer")



microATAC_trk <- DataTrack(data = MicroATAC$value, start = MicroATAC$start,
                    end = MicroATAC$start, chromosome = MicroATAC$chr,
                    genome = 'hg19', name = "Microglia ATAC", type = 'l')


PU1H3K4_trk <- DataTrack(data = PU1H3K4$value, start = PU1H3K4$start,
                           end = PU1H3K4$start, chromosome = PU1H3K4$chr,
                           genome = 'hg19', name = "Microglia H3K4me3", type = 'l')


#################################
# 3. Import ref Genes from UCSC #
#################################
#* Import Reference genes from UCSC website. 
#* Requires internet connection


# Set range
from <- 32300000
to <- 32750000

refGenes <- UcscTrack(genome = "hg19", chromosome = "chr6",
                      track = "knownGene", from = from, to = to,
                      trackType = "GeneRegionTrack", 
                      rstarts = "exonStarts", rends = "exonEnds", 
                      gene = "name",  symbol = "name2", 
                      transcript = "name", strand = "strand",
                      fill = "#8282d2", stacking = "dense", 
                      name = "Ref Genes")




#########################################
# 4. Create Chromatin interaction track #
#########################################

# Import dataset (go to URL in description for source)
microgliaPlac = read.csv('MicrogliaPlacseq.csv')

# Create Granges for Start Anchor
microgliaPlacstart = GRanges(seqnames = microgliaPlac$X.chrom,
                             ranges = IRanges(start = microgliaPlac$region1Start,
                                              end = microgliaPlac$regionEnd))

# Create Granges for end Anchor
microgliaPlacEnd = GRanges(seqnames = microgliaPlac$X.chrom,
                           ranges = IRanges(start = microgliaPlac$region2Start,
                                            end = microgliaPlac$region2End))

# Combine anchors into single Grange genomic interaction object
microgliaPlacinteractionObj <- GenomicInteractions(
  anchor1 = microgliaPlacstart,
  anchor2 = microgliaPlacEnd,
);

# Create into Gviz track
microgliaPlacInteractionTrk <- InteractionTrack(
  microgliaPlacinteractionObj, 
  name = 'Microglia PlacSeq'
);

# Ensure its Chromosome 6
microgliaPlacInteractionTrk@chromosome = 'chr6'


###################
# 5. Create Plots #
###################

# Full region all tracks -------------------------------------------------

# Get Chromosome region
ideom = IdeogramTrack(genome = 'hg19', chromosome = 'chr6')
genome.axis.track <- GenomeAxisTrack()

# create highlight track
#* overlay fine mapped variant locations
ht <- HighlightTrack(trackList = list(refGenes, LHx2En_Trk , LHx2prom_Trk, LHx2SupEn_Trk,
                                      PU1En_Trk , PU1prom_Trk, PU1SupEn_Trk,
                                      neuonEn_Trk , neuonprom_Trk, neuonSupEn_Trk,
                                      olig2En_Trk , olig2prom_Trk, olig2SupEn_Trk,
                                      microATAC_trk, PU1H3K4_trk, microgliaPlacInteractionTrk),
                     # fine-map locations
                     start = c(32431147,32441199,32544901,32550331,32550828,32550860,32560631,32561331,32572113
                               ,32572302,32572305,32575513,32576592,32576688,32577222,32578772,32583027,32583357
                               ,32584693,32606941,32606949,32607141,32607729,32610825,32611590,32619144),
                     width = 1, chromosome = 6, col = 'black')

# Plot

plotTracks(list(ideom, genome.axis.track, ht),
           from = 32300000, to = 32750000, 
           sizes = c(1,5,2,3.5,3.5,3.5,3.5,3.5,3.5,3.5,3.5,3.5,3.5,3.5,3.5,3.5,3.5,5), cex = .7)

# Zoomed in region limited tracks tracks -------------------------------------------------

# Create Highlight track with overlayed fine mapped variants
#* Select relevanr gene tracks only
ht2 <- HighlightTrack(trackList = list(refGenes,
                                      PU1En_Trk , PU1prom_Trk, PU1SupEn_Trk,
                                      microATAC_trk, PU1H3K4_trk, microgliaPlacInteractionTrk),
                     start = c(32431147,32441199,32544901,32550331,32550828,32550860,32560631,32561331,32572113
                               ,32572302,32572305,32575513,32576592,32576688,32577222,32578772,32583027,32583357
                               ,32584693,32606941,32606949,32607141,32607729,32610825,32611590,32619144),
                     width = 1, chromosome = 6, col = 'black')

# Plot with Zoomed region
plotTracks(list(ideom, genome.axis.track, ht2),
           from = 3250000, to = 32590000,
           sizes = c(1,1.2,1,2,2,2,2.2,2.2,2), cex = .8)









ht3 <- HighlightTrack(trackList = list(refGenes,
                                       PU1En_Trk , PU1prom_Trk, PU1SupEn_Trk,
                                       microATAC_trk, PU1H3K4_trk, microgliaPlacInteractionTrk),
                      start = c(32431147,32441199,32544901),
                      width = 1, chromosome = 6, col = 'black',)


plotTracks(list(ideom, genome.axis.track, ht3),
          from = 32400000, to = 32500000)





