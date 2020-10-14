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

##data frame to reduce loading
butternut_latlon <- read.csv("DataFiles/24Populations/11loci/butternut_24pops_lonlat.csv")

##first reduce it by missing data 
butternut_latlon_nomd <- butternut_latlon[butternut_latlon$Ind %in% rownames(butternutgen_nomd@tab),]

##8 loci lon lat
butternut_8loci <- read.csv("DataFiles\\24Populations\\8loci\\lon_lat_8loci.csv")

##reduce for FST file 
butternut_8loci <- butternut_8loci[,-c(1:2,4:5)]

##calculate PWFst 
fst_11loci <- as.matrix(pairwise.neifst(butternut_24pop[,-1]))

fst_8loci <- as.matrix(pairwise.neifst(butternut_8loci[,-1]))

##load in documents
fstlat_24pop_11loci <- read.csv("DataFiles\\24Populations\\11loci\\fst_lat_11loci.csv")
fstlat_24pop_8loci <- read.csv("DataFiles\\24Populations\\8loci\\fst_lat_8loci.csv")

##read in mean coordinate df
mean_coord_df <- read.csv("DataFiles\\24Populations\\mean_coords.csv")

##remove first column
mean_coord_df <- mean_coord_df[,-c(1,4)]

##create a length vector
butternut_pops_n <- length(fstlat_24pop_8loci$X)

##calculate distance to identify if there is isolation by distance 
butternut_dist <- matrix(nrow = butternut_pops_n, ncol = butternut_pops_n)

for(first in 1:butternut_pops_n){
  for(second in 1:butternut_pops_n){
    butternut_dist[first,second] <-  distm(mean_coord_df[first,], mean_coord_df[second,], fun = distGeo)/1000
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
