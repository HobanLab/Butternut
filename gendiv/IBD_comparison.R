######################
##### Libraries ######
######################

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

######################################
###### Genetic Analyses - IBD ########
######################################
butternutgen_24pop <- read.genepop("DataFiles/24Populations/11loci/butternut_24pops.gen", ncode = 3)

##reduce based on missing data 
butternutgen_nomd <- missingno(butternut_24pop_gen, type = "geno", cutoff = 0.25, quiet = FALSE, freq = FALSE)

##load in documents
butternut_24pop_fstlat <- read.csv("DataFiles\\24Populations\\11loci\\fst_lat_11loci.csv")


##create a length vector
butternut_pops_n <- length(unique(butternut_24pop$Pop))

##calculate distance to identify if there is isolation by distance 
butternut_dist <- matrix(nrow = butternut_pops_n, ncol = butternut_pops_n)

for(first in 1:butternut_pops_n){
  for(second in 1:butternut_pops_n){
    butternut_dist[first,second] <-  distm(butternut_coord_df[first,], butternut_coord_df[second,], fun = distGeo)/1000
  }
}

##generate IBD plot

##11 loci IBD plot
dist_fst_11loci_24pop <- array(c(butternut_24pop_fst,butternut_dist), dim = c(butternut_pops_n,butternut_pops_n,2))
pdf("dist_fst_11loci_24pop.pdf", width=10,height=8)
plot(dist_fst_11loci_24pop[12:24,1:11,2],dist_fst_11loci_24pop[12:24,1:11,1], pch = 16, col = "mediumpurple1", xlab = "Distance (km)", ylab = "PWFst")
points(dist_fst_11loci_24pop[1:11,1:11,2],dist_fst_11loci_24pop[1:11,1:11,1], col = "red1", xlab = "Distance (km)", ylab = "PWFst", pch = 16)
points(dist_fst_11loci_24pop[12:24,12:24,2],dist_fst_11loci_24pop[12:24,12:24,1], pch = 16, col = "dodgerblue4")
legend('topright',c("Canada vs Canada", "US v US", "US v Canada"), pch = 16, col = c("red1","dodgerblue4","mediumpurple1"), cex = 0.75, pt.cex = 0.75)
dev.off()