##########################
######## Libraries #######
##########################

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

###############################################
############## Working Directory ##############
###############################################

setwd("G:/My Drive/Hoban_Lab_Docs/Projects/Butternut_JUCI")

###############################################
############# Individual Reduction ############
###############################################

##following pop reclassification - load 24 pop document
butternut_24pop_gen <- read.genepop("G:/My Drive/Hoban_Lab_Docs/Projects/Butternut_JUCI/DataFiles/24Populations/11loci/butternut_24pops.gen", ncode = 3)

##now remove individuals with greater than 25% missing data
butternutgen_nomd <- missingno(butternut_24pop_gen, type = "geno", cutoff = 0.25, quiet = FALSE, freq = FALSE)

#load in relatedness file
butternut_24pop_rel <- read.csv("G:/My Drive/Hoban_Lab_Docs/Projects/Butternut_JUCI/DataFiles/24Populations/11loci/butternut_24pop_rel.csv")

##calculate relatedness
butternutrel <- Demerelate(butternut_24pop_rel, object = T, value = "loiselle")

##now identify how many individuals have greater than 25% relatedness = half siblings
butternut_halfsib_names <- names(which(unlist(butternutrel$Empirical_Relatedness) > 0.25))

##then use this to create a document which has all the unique individual numbers for every highly related individuals
butternut_halfsib_names_cleanfront <- gsub("^.*\\.","", butternut_halfsib_names)

butternut_halfsib_names_cleanback <- gsub("^.*\\_","", butternut_halfsib_names_cleanfront)

relate_ind_remove <- unique(butternut_halfsib_names_cleanback)

##then subset genind file
butternutgen_24pop_reduced <- butternutgen_nomd[!rownames(butternutgen_nomd@tab) %in% relate_ind_remove,]

##data frame to reduce loading
butternut_latlon <- read.csv("G:/My Drive/Hoban_Lab_Docs/Projects/Butternut_JUCI/DataFiles/24Populations/11loci/butternut_24pops_lonlat.csv")

##first reduce it by missing data 
butternut_latlon_nomd <- butternut_latlon[butternut_latlon$Ind %in% butternut_24pop_rel$Ind,]

##then reduce it by relatedness individuals
butternut_24pop_coords <- butternut_latlon_nomd[!butternut_latlon_nomd$Ind %in% relate_ind_remove,]