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
library(dplyr)

###############################################
############## Working Directory ##############
###############################################

setwd("G:/My Drive/Hoban_Lab_Docs/Projects/Butternut_JUCI")

######################################
###### Genetic Analyses - IBD ########
######################################
##data frame for fst
butternut_latlon <- read.csv("DataFiles/24Populations/11loci/butternut_24pops_lonlat.csv")

##prepare fst doc 
fst_11loci_ready <- butternut_latlon[,-c(1,3:4)]

##first reduce it by missing data 
butternut_latlon_nomd <- butternut_latlon[butternut_latlon$Ind %in% rownames(butternutgen_nomd@tab),]

##8 loci lon lat
butternut_8loci <- read.csv("DataFiles\\24Populations\\8loci\\butternut_24pop_loc.csv")

##Prep 8 loci doc 
fst_8loci_ready <- butternut_8loci[,-1]

##calculate PWFst 
fst_11loci <- as.matrix(pairwise.neifst(fst_11loci_ready))

fst_8loci <- as.matrix(pairwise.neifst(fst_8loci_ready))


##read in mean coordinate df
mean_coord_df <- read.csv("DataFiles\\24Populations\\mean_coords.csv")

##reorder to match Fst
rownames(mean_coord_df) <- mean_coord_df[,1]
mean_coord_df <- mean_coord_df[rownames(fst_11loci),]

##remove first column
mean_coord_df <- mean_coord_df[,-1]

##create a length vector
butternut_pops_n <- length(mean_coord_df$Mean.Lon)

##calculate distance to identify if there is isolation by distance 
butternut_dist <- matrix(nrow = butternut_pops_n, ncol = butternut_pops_n)

for(first in 1:butternut_pops_n){
  for(second in 1:butternut_pops_n){
    butternut_dist[first,second] <-  distm(mean_coord_df[first,], mean_coord_df[second,], fun = distGeo)/1000
  }
}

#######generate IBD plot

##11 loci IBD plot
IBD_11loci_array <- array(c(fst_11loci,butternut_dist), dim = c(butternut_pops_n,butternut_pops_n,2))

##create linear regression
IBD_11loci_lm <- lm(as.numeric(IBD_11loci_array[,,1])~as.numeric(IBD_11loci_array[,,2]))
IBD_11loci_lm_summary <- summary(IBD_11loci_lm)

##r2 creation 
r2_IBD <- IBD_11loci_lm_summary$adj.r.squared
coef_IBD <- IBD_11loci_lm_summary$coefficients
pvalue_IBD <- coef_IBD[2,4]
rp_IBD <- vector('expression',2)
rp_IBD[1] <- substitute(expression(italic(R)^2 == MYVALUE), 
                            list(MYVALUE = format(r2_IBD,dig=3)))[2]
rp_IBD[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                           list(MYOTHERVALUE = format(pvalue_IBD, digits = 2)))[2]

##create pdf
pdf("Graphical_Stat_Results\\PostIndRemoval\\24pop\\11loci\\gendiv\\IBD_11loci_24pop.pdf", width=10,height=8)

##plot docs
plot(IBD_11loci_array[12:24,1:11,2],IBD_11loci_array[12:24,1:11,1], pch = 16, col = "mediumpurple1", xlab = "Distance (km)", ylab = "PWFst", main = "11 loci IBD")
points(IBD_11loci_array[1:11,1:11,2],IBD_11loci_array[1:11,1:11,1], col = "red1", xlab = "Distance (km)", ylab = "PWFst", pch = 16)
points(IBD_11loci_array[12:24,12:24,2],IBD_11loci_array[12:24,12:24,1], pch = 16, col = "dodgerblue4")
legend('topright',c("Canada vs Canada", "US v US", "US v Canada"), pch = 16, col = c("red1","dodgerblue4","mediumpurple1"), cex = 0.75, pt.cex = 0.75)
legend('topleft', legend = rp_IBD, bty = 'n', border = "black", pt.cex = 1, cex = 0.8)
abline(IBD_11loci_lm, col = "blue")
dev.off()

###8 loci IBD plot
IBD_8loci_array <- array(c(fst_8loci,butternut_dist), dim = c(butternut_pops_n,butternut_pops_n,2))

##create linear regression
IBD_8loci_lm <- lm(as.numeric(IBD_8loci_array[,,1])~as.numeric(IBD_8loci_array[,,2]))
IBD_8loci_lm_summary <- summary(IBD_8loci_lm)

##r2 creation 
r2_8loci_IBD <- IBD_8loci_lm_summary$adj.r.squared
coef_8loci_IBD <- IBD_8loci_lm_summary$coefficients
pvalue_8loci_IBD <- coef_8loci_IBD[2,4]
rp_8loci_IBD <- vector('expression',2)
rp_8loci_IBD[1] <- substitute(expression(italic(R)^2 == MYVALUE), 
                        list(MYVALUE = format(r2_8loci_IBD,dig=3)))[2]
rp_8loci_IBD[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                       list(MYOTHERVALUE = format(pvalue_8loci_IBD, digits = 2)))[2]

##create pdf
pdf("Graphical_Stat_Results\\PostIndRemoval\\24pop\\8loci\\gendiv\\IBD_8loci_24pop.pdf", width=10,height=8)

##plot docs
plot(IBD_8loci_array[12:24,1:11,2],IBD_8loci_array[12:24,1:11,1], pch = 16, col = "mediumpurple1", xlab = "Distance (km)", ylab = "PWFst", main = "IBD 8 loci", ylim = c(0,0.07))
points(IBD_8loci_array[1:11,1:11,2],IBD_8loci_array[1:11,1:11,1], col = "red1", xlab = "Distance (km)", ylab = "PWFst", pch = 16)
points(IBD_8loci_array[12:24,12:24,2],IBD_8loci_array[12:24,12:24,1], pch = 16, col = "dodgerblue4")
legend('topright',c("Canada vs Canada", "US v US", "US v Canada"), pch = 16, col = c("red1","dodgerblue4","mediumpurple1"), cex = 0.75, pt.cex = 0.75)
legend('topleft', legend = rp_8loci_IBD, bty = 'n', border = "black", pt.cex = 1, cex = 0.8)
abline(IBD_8loci_lm, col = "blue")
dev.off()

