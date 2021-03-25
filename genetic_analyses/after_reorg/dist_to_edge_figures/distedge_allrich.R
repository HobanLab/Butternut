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
dist_edge_path <- "G:\\Shared drives\\Emily_Schumacher\\Butternut\\butternut_publication_figures\\dist_edge\\"
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

######################################################################
#################### Run Analyses and Create DF  #####################
######################################################################
butternut_reorg_allrich <- colSums(allelic.richness(butternutgen_reorg)$Ar)/length(butternutgen_reorg@loc.n.all)

##create df with allelic richness and distance
all_rich_dist_df <- data.frame(dist_edge_df[,4], butternut_reorg_allrich,dist_edge_df[,5])
colnames(all_rich_dist_df) <- c("Distance", "Allelic_Richness", "Col")
rownames(all_rich_dist_df) <- butternut_24pop_names

##reduced data frame
allrich_red <- all_rich_dist_df[-c(16,22,23),]

######################################################################
########################## Linear Models  ############################
######################################################################

###now calculate regressions - linear
allrich_lm <- lm(all_rich_dist_df[,2]~all_rich_dist_df[,1])
allrich_lm_sum <- summary(allrich_lm)

##Reduced data frame
allrich_red_lm <- lm(allrich_red[,2]~allrich_red[,1])
allrich_red_lm_sum <- summary(allrich_red_lm)

##create r2 and p for full df 
linear_allrich_dist_r2 <- allrich_lm_sum$adj.r.squared
linear_allrich_dist_coef <- allrich_lm_sum$coefficients
linear_allrich_dist_pvalue <- linear_allrich_dist_coef[2,4]
linear_allrich_dist_rp <- vector('expression',2)
linear_allrich_dist_rp[1] <- substitute(expression(italic(R)^2 == MYVALUE), 
                            list(MYVALUE = format(linear_allrich_dist_r2,dig=3)))[2]
linear_allrich_dist_rp[2] <- substitute(expression(italic(p) == MYOTHERVALUE), 
                            list(MYOTHERVALUE = format(linear_allrich_dist_pvalue, digits = 2)))[2]

##create r2 and p for reduced df 
linear_allrich_dist_red_r2 <- allrich_red_lm_sum$adj.r.squared
linear_allrich_dist_red_coef <- allrich_red_lm_sum$coefficients
linear_allrich_dist_red_pvalue <- linear_allrich_dist_red_coef[2,4]
linear_allrich_dist_red_rp <- vector('expression',2)
linear_allrich_dist_red_rp[1] <- substitute(expression(italic(R)^2 == MYVALUE), 
                                        list(MYVALUE = format(linear_allrich_dist_red_r2,dig=3)))[2]
linear_allrich_dist_red_rp[2] <- substitute(expression(italic(p) == MYOTHERVALUE), 
                                        list(MYOTHERVALUE = format(linear_allrich_dist_red_pvalue, digits = 2)))[2]

##################################################################
################# Plotting Linear Models  ########################
##################################################################

##Allelic Richness Distance
pdf(paste0(dist_edge_path, "linear_allrich_dist.pdf"), width = 8, height = 6)

plot(all_rich_dist_df[,2]~all_rich_dist_df[,1], 
     col = as.character(all_rich_dist_df[,3]), pch = 17,
     ylab = "Number of Alleles", xlab = "Distance to Range Edge (km)", 
     cex = (butternut_24pop_poppr[1:24,2]/42), ylim = c(5,10), xlim = c(0,600),
     main = "Number of Alleles Compared with Distance to Range Edge (km)")

text(all_rich_dist_df[,2]~all_rich_dist_df[,1], 
     labels = c(1:24), cex = 0.8, pos = 1)

abline(allrich_lm, col = "dodgerblue4", lwd = 3)
abline(allrich_red_lm, col = "darkorchid", lwd = 3)

##legends
legend('topleft', legend = linear_allrich_dist_rp, pch = 4, col = "dodgerblue4",
       title = "Regression with WI",
       bty = 'n', border = "black", pt.cex = 1, cex = 0.8)

legend('bottomleft', legend = linear_allrich_dist_red_rp, pch = 4, 
       col = "darkorchid", title = "Regression without WI",
       bty = 'n', border = "black", pt.cex = 1, cex = 0.8)

legend('bottom', legend = c("New Brunswick", "Ontario", "United States"), pch = 17, col = c("firebrick1", "firebrick4","dodgerblue"))

dev.off()

##create data frame of pvalues r2 for linear and quad models
dist_edge_allrich_df <- matrix(nrow = 4, ncol = 2)
rownames(dist_edge_allrich_df) <- c("linear_allrich_dist_edge",
                                    "linear_allrich_dist_edge_woWI",
                                    "quadratic_allrich_dist_edge",
                                    "quadratic_allrich_dist_edge_woWI")
colnames(dist_edge_allrich_df) <- c("Pvalue","R2")
##write values
dist_edge_allrich_df[1,] <- c(linear_allrich_dist_pvalue, linear_allrich_dist_r2)
dist_edge_allrich_df[2,] <- c(linear_allrich_dist_red_pvalue, linear_allrich_dist_red_r2)
dist_edge_allrich_df[3,] <- c(quad_allrich_dist_pvalue, quad_allrich_dist_r2)
dist_edge_allrich_df[4,] <- c(quad_allrich_dist_red_pvalue, quad_allrich_dist_red_r2)

##write out
write.csv(dist_edge_allrich_df, paste0("G:\\Shared drives\\Emily_Schumacher\\butternut_publication_figures\\pvalue_r2_dist_edge_allrich.csv"))

