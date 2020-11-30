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

#####################################
############ Directories ############
#####################################
dist_edge_path <- "G:\\Shared drives\\Emily_Schumacher\\butternut_publication_figures\\dist_edge\\"
butternut_drive <- "G:\\My Drive\\Hoban_Lab_Docs\\Projects\\Butternut_JUCI\\"

##load butternut
setwd(dist_edge_path)
butternut_buffer <- shapefile("butternut_buffer")

##############################################################
#################### Load Genetic Files  #####################
##############################################################
setwd(butternut_drive)

##load current working reorganized genepop file 
butternutgen_reorg <- read.genepop("DataFiles\\24Populations\\reorg\\reorg_gen_24pop.gen", ncode = 3)

##Lon lat files
butternut_reorg_lonlat <- read.csv("DataFiles\\24Populations\\reorg\\reorg_lon_lat.csv")

###Name the reorg file
rownames(butternutgen_reorg@tab) <- butternut_reorg_lonlat$Ind

##GESTE Results
reorg_geste_fst <- read.csv(paste0(butternut_drive,"\\DataFiles\\24Populations\\reorg\\GESTE_fst.csv"))

##create population name file 
butternut_24pop_names <- unique(butternut_reorg_lonlat$Pop)

##Name levels
levels(butternutgen_reorg@pop) <- butternut_24pop_names

##calculate poppr
butternut_24pop_poppr <- poppr(butternutgen_reorg)

##get distance to range edge for each pop
dist_edge_df <- read.csv(paste0(dist_edge_path, "dist_edge_df.csv"))

#############################################################################################
#################### Run Analyses and Create DF of Just Ontario and US  #####################
#############################################################################################
butternut_reorg_allrich <- colSums(allelic.richness(butternutgen_reorg)$Ar)/length(butternutgen_reorg@loc.n.all)

##create df with allelic richness and distance
all_rich_dist_df <- data.frame(dist_edge_df[,4], butternut_reorg_allrich,dist_edge_df[,5])
colnames(all_rich_dist_df) <- c("Distance", "Allelic_Richness", "Col")
rownames(all_rich_dist_df) <- butternut_24pop_names

##remove New Brunswick
all_rich_dist_df_wo_NB <- all_rich_dist_df[-c(1:6),]
all_rich_dist_df_wo_NB_red <- all_rich_dist_df_wo_NB[-c(10,16:17),]

######################################################################
########################## Linear Models  ############################
######################################################################

###now calculate regressions - linear
allrich_wo_NB_lm <- lm(all_rich_dist_df_wo_NB[,2]~all_rich_dist_df_wo_NB[,1])
allrich_wo_NB_lm_sum <- summary(allrich_wo_NB_lm)

##Reduced data frame
allrich_red_lm <- lm(all_rich_dist_df_wo_NB_red[,2]~all_rich_dist_df_wo_NB_red[,1])
allrich_red_lm_sum <- summary(allrich_red_lm)

##create r2 and p for full df 
linear_allrich_dist_wo_NB_r2 <- allrich_wo_NB_lm_sum$adj.r.squared
linear_allrich_dist_wo_NB_coef <- allrich_wo_NB_lm_sum$coefficients
linear_allrich_dist_wo_NB_pvalue <- linear_allrich_dist_wo_NB_coef[2,4]
linear_allrich_dist_wo_NB_rp <- vector('expression',2)
linear_allrich_dist_wo_NB_rp[1] <- substitute(expression(italic(R)^2 == MYVALUE), 
                                        list(MYVALUE = format(linear_allrich_dist_wo_NB_r2,dig=3)))[2]
linear_allrich_dist_wo_NB_rp[2] <- substitute(expression(italic(p) == MYOTHERVALUE), 
                                        list(MYOTHERVALUE = format(linear_allrich_dist_wo_NB_pvalue, digits = 2)))[2]

##create r2 and p for reduced df 
linear_allrich_dist_wo_NB_red_r2 <- allrich_red_lm_sum$adj.r.squared
linear_allrich_dist_wo_NB_red_coef <- allrich_red_lm_sum$coefficients
linear_allrich_dist_wo_NB_red_pvalue <- linear_allrich_dist_wo_NB_red_coef[2,4]
linear_allrich_dist_wo_NB_red_rp <- vector('expression',2)
linear_allrich_dist_wo_NB_red_rp[1] <- substitute(expression(italic(R)^2 == MYVALUE), 
                                          list(MYVALUE = format(linear_allrich_dist_wo_NB_red_r2,dig=3)))[2]
linear_allrich_dist_wo_NB_red_rp[2] <- substitute(expression(italic(p) == MYOTHERVALUE), 
                                         list(MYOTHERVALUE = format(linear_allrich_dist_wo_NB_red_pvalue, digits = 2)))[2]

##################################################################
################# Plotting Linear Models  ########################
##################################################################

##Allelic Richness Distance
pdf(paste0(dist_edge_path, "linear_allrich_dist_wo_NB.pdf"), width = 8, height = 6)

plot(all_rich_dist_df_wo_NB[,2]~all_rich_dist_df_wo_NB[,1], 
     col = all_rich_dist_df_wo_NB[,3], pch = 17,
     ylab = "Number of Alleles", xlab = "Distance to Range Edge (km)", 
     cex = (butternut_24pop_poppr[7:24,2]/80), ylim = c(5,10), xlim = c(0,600),
     main = "Number of Alleles Compared with Distance to Range Edge (km)")

text(all_rich_dist_df_wo_NB[,2]~all_rich_dist_df_wo_NB[,1], 
     labels = rownames(all_rich_dist_df_wo_NB), cex = 0.8, pos = 1)

abline(allrich_wo_NB_lm, col = "dodgerblue4")
abline(allrich_red_lm, col = "darkorchid")

##legends
legend('topleft', legend = linear_allrich_dist_wo_NB_rp, pch = 4, col = "dodgerblue4",
       title = "Regression with WI",
       bty = 'n', border = "black", pt.cex = 1, cex = 0.8)

legend('bottomleft', legend = linear_allrich_dist_wo_NB_red_rp, pch = 4, 
     col = "darkorchid", title = "Regression without WI",
     bty = 'n', border = "black", pt.cex = 1, cex = 0.8)

legend('bottom', legend = c("Ontario", "United States"), pch = 17, col = c("firebrick4","dodgerblue"))

dev.off()

###########################################################
################### Create DF of Just NB ##################
###########################################################
##create df with allelic richness and distance
NB_all_rich_dist_df <- all_rich_dist_df[1:6,]

##################################################################
################# Plotting Linear Models  ########################
##################################################################
###now calculate regressions - linear
allrich_NB_lm <- lm(NB_all_rich_dist_df[,2]~NB_all_rich_dist_df[,1])
allrich_NB_lm_sum <- summary(allrich_NB_lm)

##create r2 and p for full df 
linear_allrich_NB_r2 <- allrich_NB_lm_sum$adj.r.squared
linear_allrich_NB_coef <- allrich_NB_lm_sum$coefficients
linear_allrich_NB_pvalue <- linear_allrich_NB_coef[2,4]
linear_allrich_NB_rp <- vector('expression',2)
linear_allrich_NB_rp[1] <- substitute(expression(italic(R)^2 == MYVALUE), 
                                              list(MYVALUE = format(linear_allrich_NB_r2,dig=3)))[2]
linear_allrich_NB_rp[2] <- substitute(expression(italic(p) == MYOTHERVALUE), 
                                              list(MYOTHERVALUE = format(linear_allrich_NB_pvalue, digits = 2)))[2]


##################################################################
################# Plotting Linear Models  ########################
##################################################################

##Allelic Richness Distance
pdf(paste0(dist_edge_path, "linear_allrich_dist_NB.pdf"), width = 8, height = 6)

plot(NB_all_rich_dist_df[,2]~NB_all_rich_dist_df[,1], 
     col = NB_all_rich_dist_df[,3], pch = 17,
     ylab = "Number of Alleles", xlab = "Distance to Range Edge (km)", 
     cex = (butternut_24pop_poppr[7:24,2]/80), ylim = c(5,10), xlim = c(0,300),
     main = "Number of Alleles Compared with Distance to Range Edge (km)")

text(NB_all_rich_dist_df[,2]~NB_all_rich_dist_df[,1], 
     labels = rownames(NB_all_rich_dist_df), cex = 0.8, pos = 1)

abline(allrich_NB_lm, col = "dodgerblue4")
#abline(allrich_red_lm, col = "darkorchid")

##legends
legend('topleft', legend = linear_allrich_NB_rp, pch = 4, col = "dodgerblue4",
       title = "New Brunswick",
       bty = 'n', border = "black", pt.cex = 1, cex = 0.8)

#legend('bottomleft', legend = linear_allrich_dist_wo_NB_red_rp, pch = 4, 
 #      col = "darkorchid", title = "Regression without WI",
   #    bty = 'n', border = "black", pt.cex = 1, cex = 0.8)

legend('bottom', legend = c("New Brunswick"), pch = 17, col = c("firebrick1"))

dev.off()
