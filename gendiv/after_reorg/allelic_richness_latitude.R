##########################
######## Libraries #######
##########################

library(diveRsity)
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
library(plotrix)
library(ggpmisc)
library(factoextra)
library(GISTools)
library(raster)
library(rgdal)
library(sp)
library(PopGenReport)

#####################################
############ Directories ############
#####################################
shared_drive <- "G:\\Shared drives\\Emily_Schumacher\\butternut_publication_figures"
butternut_drive <- "G:\\My Drive\\Hoban_Lab_Docs\\Projects\\Butternut_JUCI"

#####################################
############# Load Files ############
#####################################
setwd(paste0(butternut_drive))

##load current working reorganized genepop file 
butternutgen_reorg <- read.genepop("DataFiles\\24Populations\\reorg\\reorg_gen_24pop.gen", ncode = 3)

##load relatedness
reorg_relatedness <- read.csv("DataFiles\\24Populations\\reorg\\reorg_relatedness.csv")

##rename genind 
rownames(butternutgen_reorg@tab) <- reorg_relatedness$Ind

##Lon lat files 
butternut_reorg_lonlat <- read.csv("DataFiles\\24Populations\\reorg\\reorg_lon_lat.csv")

##pop name
butternut_24pop_names <- unique(butternut_reorg_lonlat$Pop)

##create butternut poppr
butternut_poppr <- poppr(butternutgen_reorg)

#############################################
############# Mean Lat and Long #############
#############################################

##Calculate mean latitude and longitude for every population
##now do lon lat calcs
butternut_mean_lon <- matrix()
butternut_mean_lat <- matrix()

##loops for mean lat/lon
for(pop in butternut_24pop_names){
  
  butternut_mean_lon[pop] <- mean(butternut_reorg_lonlat[butternut_reorg_lonlat$Pop == pop,][,3])  
  
}

for(pop in butternut_24pop_names){
  
  butternut_mean_lat[pop] <- mean(butternut_reorg_lonlat[butternut_reorg_lonlat$Pop == pop,][,4])  
  
}

##convert to matrix
butternut_mean_lon <- matrix(butternut_mean_lon)
butternut_mean_lat <- matrix(butternut_mean_lat)

##document cleanup
butternut_mean_lon <- butternut_mean_lon[-1]
butternut_mean_lat <- butternut_mean_lat[-1]

##create min and max documents
min_lon <- min(butternut_mean_lon)
max_lon <- max(butternut_mean_lon)
min_lat <- min(butternut_mean_lat)
max_lat <- max(butternut_mean_lat)

##########################################################
################ Genetic Diversity Plotting ##############
##########################################################
################Allelic Richness
##Calculate Allelic Richness
reorg_allrich <- colSums(allelic.richness(butternutgen_reorg)$Ar)/length(butternutgen_reorg@loc.n.all)

##create df with allelic richness and latitude
all_rich_lat_df <- data.frame(butternut_mean_lat, reorg_allrich)
colnames(all_rich_lat_df) <- c("Mean_Lat", "Allelic_Richness")
rownames(all_rich_lat_df) <- butternut_24pop_names
all_rich_lat_df$Color <- NA
all_rich_lat_df[1:6,3] <- "firebrick1"
all_rich_lat_df[7:11,3] <- "firebrick4"
all_rich_lat_df[12:24,3] <- "dodgerblue"

##Create linear relationship with lat and number of alleles 
allrich_lm <- lm(all_rich_lat_df[,2]~all_rich_lat_df[,1])
allrich_lm_sum <- summary(allrich_lm)

###df reduced without wisconsin
allrich_lat_red_df <- all_rich_lat_df[-c(16,22,23),]

##Now create linear model 
allrich_red_lm <- lm(allrich_lat_red_df[,2]~allrich_lat_red_df[,1])
allrich_red_lm_sum <- summary(allrich_red_lm)

##create r2 and p values
allrich_r2 <- allrich_lm_sum$adj.r.squared
allrich_coef <- allrich_lm_sum$coefficients
allrich_pvalue <- allrich_coef[2,4]
allrich_rp <- vector('expression',2)
allrich_rp[1] <- substitute(expression(italic(R)^2 == MYVALUE), 
                                list(MYVALUE = format(allrich_r2,dig=3)))[2]
allrich_rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                               list(MYOTHERVALUE = format(allrich_pvalue, digits = 2)))[2]

##create r2 and p values for reduced df 
allrich_red_r2 <- allrich_red_lm_sum$adj.r.squared
allrich_red_coef <- allrich_red_lm_sum$coefficients
allrich_red_pvalue <- allrich_red_coef[2,4]
allrich_red_rp <- vector('expression',2)
allrich_red_rp[1] <- substitute(expression(italic(R)^2 == MYVALUE), 
                            list(MYVALUE = format(allrich_red_r2,dig=3)))[2]
allrich_red_rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                           list(MYOTHERVALUE = format(allrich_red_pvalue, digits = 2)))[2]

##Allelic Richness Linear
setwd("G:\\Shared drives\\Emily_Schumacher\\butternut_publication_figures")
pdf("AllRich_Lat.pdf", width = 8, height = 6)
plot(all_rich_lat_df[,2]~all_rich_lat_df[,1], col = all_rich_lat_df[,3], pch = 17, main = "Number of Alleles Compared to Mean Latitude", ylab = "Number of Alleles", xlab = "Mean Latitude", cex = (butternut_poppr[1:24,2]/80), ylim = c(5,10))
text(all_rich_lat_df[,2]~all_rich_lat_df[,1], labels = butternut_24pop_names, cex = 0.8, pos = 1)
abline(allrich_lm, col = "dodgerblue4")
abline(allrich_red_lm, col = "darkorchid4")
legend('topleft', legend = allrich_rp, bty = 'n', border = "black", pt.cex = 1, cex = 1, pch = 4, col = "dodgerblue4")
legend('bottomleft', legend = allrich_red_rp, bty = 'n', border = "black", pt.cex = 1, cex = 1, pch = 4, col = "darkorchid4")
legend('bottom', legend = c("New Brunswick", "Ontario", "United States"), pch = 17, col = c("firebrick1", "firebrick4","dodgerblue"))
dev.off()

##Do quadratic
#calculate data points 
points <- all_rich_lat_df[,1]
points2 <- points^2

###Quadratic linear model
quadratic_model_lm <-lm(all_rich_lat_df[,2] ~ 
                          points + points2)
##Calculate summary documents 
quadratic_model_lm_sum <- summary(quadratic_model_lm)

#calculate model out of points 
points_values <- seq(min_lat, max_lat, 0.5)
points_counts <- predict(quadratic_model_lm,list(points=points_values, points2=points_values^2))

###pointsance without wisconsin populations
points_21 <- allrich_lat_red_df[,1]
points2_21 <- allrich_lat_red_df[,1]^2

##create quadratic regression documents 
quadratic_21pops_model_lm <-lm(allrich_lat_red_df[,2] ~ 
                                 points_21 + points2_21)
quadratic_21pops_model_lm_sum <- summary(quadratic_21pops_model_lm)

##values and pointsance
points_values_21 <- seq(min_lat, max_lat, 0.5)
points_counts_21 <- predict(quadratic_21pops_model_lm,list(points_21=points_values_21, points2_21=points_values_21^2))

##create r2 and p values
allrich_quad_r2 <- quadratic_model_lm_sum$adj.r.squared
allrich_quad_coef <- quadratic_model_lm_sum$coefficients
allrich_quad_pvalue <- allrich_quad_coef[2,4]
allrich_quad_rp <- vector('expression',2)
allrich_quad_rp[1] <- substitute(expression(italic(R)^2 == MYVALUE), 
                            list(MYVALUE = format(allrich_quad_r2,dig=3)))[2]
allrich_quad_rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                           list(MYOTHERVALUE = format(allrich_quad_pvalue, digits = 2)))[2]

##create r2 and p values for reduced df 
allrich_quad_red_r2 <- quadratic_21pops_model_lm_sum$adj.r.squared
allrich_quad_red_coef <- quadratic_21pops_model_lm_sum$coefficients
allrich_quad_red_pvalue <- allrich_quad_red_coef[2,4]
allrich_quad_red_rp <- vector('expression',2)
allrich_quad_red_rp[1] <- substitute(expression(italic(R)^2 == MYVALUE), 
                                list(MYVALUE = format(allrich_quad_red_r2,dig=3)))[2]
allrich_quad_red_rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                               list(MYOTHERVALUE = format(allrich_quad_red_pvalue, digits = 2)))[2]

##Now Plot 
pdf("quad_reg_num_all.pdf", width = 8, height = 6)
plot(all_rich_lat_df[,2]~all_rich_lat_df[,1], col = all_rich_lat_df[,3], 
     pch = 17, main = "Number of Alleles Compared to Mean Latitude", ylab = "Number of Alleles", 
     xlab = "Mean Latitude", 
     cex = (butternut_poppr[1:24,2]/100), ylim = c(5,10))
text(all_rich_lat_df[,2]~all_rich_lat_df[,1], labels = butternut_24pop_names, cex = 0.8, pos = 1)
lines(points_values, points_counts, col = "darkslategray3", lwd = 3)
lines(points_values_21, points_counts_21, col = "darkseagreen4", lwd = 3)
legend('bottom', legend = c("New Brunswick", "Ontario", "United States"), pch = 17, 
       col = c("firebrick1", "firebrick4","dodgerblue"))
legend('topleft', legend = allrich_quad_rp, bty = 'n', border = "black", pt.cex = 1, cex = 1, pch = 4, col = "darkslategray3")
legend('bottomleft', legend = allrich_quad_red_rp, bty = 'n', border = "black", pt.cex = 1, cex = 1, pch = 4, col = "darkseagreen4")
dev.off()

