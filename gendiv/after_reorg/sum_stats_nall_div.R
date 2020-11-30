######################################
########### Load Libraries ###########
######################################

library(adegenet)
library(stringr)
library(tidyr)
library(hierfstat)
library(poppr)
library(Demerelate)
library(rworldmap)
library(data.table)
library(ggplot2)
library(ggrepel)
library(geosphere)
library(strataG)
library(plotrix)
library(ggpmisc)
library(factoextra)
library(diveRsity)
library(PopGenReport)

#####################################
############# Load Files ############
#####################################
##set working directory
setwd("G:/My Drive/Hoban_Lab_Docs/Projects/Butternut_JUCI")

##load reorg genind
butternutgen_reorg <- read.genepop("DataFiles/24Populations/reorg/reorg_gen_24pop.gen", ncode = 3)

##load lat lon doc 
butternut_reorg_lonlat <- read.csv("DataFiles/24Populations/reorg/reorg_lon_lat.csv")

##create name doc
butternut_24pop_names <- unique(butternut_reorg_lonlat$Pop)

##name genind doc 
rownames(butternutgen_reorg@tab) <- butternut_reorg_lonlat$Ind

######################################
############# basic stats ############
######################################
##reorg data file 
bn_sumstats <- summary(butternutgen_reorg)

BN_poppr <- poppr(butternutgen_reorg)

##
BN_hexp <- BN_poppr[1:24, 10]
bn_hexp_df <- data.frame(butternut_24pop_names,BN_hexp)


##null alleles 
bn_null_all <- null.all(butternutgen_reorg)

##create null allele table
bn_null_all_df <- matrix(nrow = length(rownames(bn_null_all$null.allele.freq$summary1)),
                         ncol = length(colnames(bn_null_all$null.allele.freq$summary1)))

bn_null_all_df <- bn_null_all$null.allele.freq$summary1

##bn HWE test
bn_hwe <- hw.test(butternutgen_reorg)

##test for LD
bn_loci <- genind2loci(butternutgen_reorg)

bn_loci <- na.omit(bn_loci)

LD2(bn_loci, details = FALSE)


