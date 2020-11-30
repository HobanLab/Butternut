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
##create hexp df 
hexp_dist_df <- data.frame(dist_edge_df[,4],as.numeric(butternut_24pop_poppr[1:24,10]), dist_edge_df[,5])
colnames(hexp_dist_df) <- c("Distance", "HExp", "Col")
rownames(hexp_dist_df) <- butternut_24pop_names

##reduced data frame
hexp_wo_wi_df <- hexp_dist_df[-c(16,22,23),]

######################################################################
########################## Linear Models  ############################
######################################################################

###now calculate regressions - linear
hexp_dist_lm <- lm(hexp_dist_df[,2]~hexp_dist_df[,1])
hexp_dist_lm_sum <- summary(hexp_dist_lm)

##linear regression without 
hexp_wo_wi_lm <- lm(hexp_wo_wi_df[,2]~hexp_wo_wi_df[,1])
hexp_wo_wi_lm_sum <- summary(hexp_wo_wi_lm)

##create r2 and p for full df 
hexp_dist_r2 <- hexp_dist_lm_sum$adj.r.squared
hexp_dist_coef <- hexp_dist_lm_sum$coefficients
hexp_dist_pvalue <- hexp_dist_coef[2,4]
hexp_dist_rp <- vector('expression',2)
hexp_dist_rp[1] <- substitute(expression(italic(R)^2 == MYVALUE), 
                                        list(MYVALUE = format(hexp_dist_r2,dig=3)))[2]
hexp_dist_rp[2] <- substitute(expression(italic(p) == MYOTHERVALUE), 
                                        list(MYOTHERVALUE = format(hexp_dist_pvalue, digits = 2)))[2]

##create r2 and p for reduced df 
hexp_wo_wi_dist_r2 <- hexp_wo_wi_lm_sum$adj.r.squared
hexp_wo_wi_dist_coef <- hexp_wo_wi_lm_sum$coefficients
hexp_wo_wi_dist_pvalue <- hexp_wo_wi_dist_coef[2,4]
hexp_wo_wi_dist_rp <- vector('expression',2)
hexp_wo_wi_dist_rp[1] <- substitute(expression(italic(R)^2 == MYVALUE), 
                                            list(MYVALUE = format(hexp_wo_wi_dist_r2,dig=3)))[2]
hexp_wo_wi_dist_rp[2] <- substitute(expression(italic(p) == MYOTHERVALUE), 
                                            list(MYOTHERVALUE = format(hexp_wo_wi_dist_pvalue, digits = 2)))[2]

##################################################################
######################## Plotting Models  ########################
##################################################################

##HExp
pdf(paste0(dist_edge_path, "hexp_dist.pdf"), width = 8, height = 6)

plot(hexp_dist_df[,2]~hexp_dist_df[,1], 
     col = hexp_dist_df[,3], pch = 17, ylim = c(0.74,0.86),
     ylab = "Expected Heterozygosity", xlab = "Distance to Range Edge (km)", 
     cex = (butternut_24pop_poppr[1:24,2]/80), xlim = c(0,600),
     main = "Expected Heterozygosity Compared with Distance to Range Edge (km)")

text(hexp_dist_df[,2]~hexp_dist_df[,1], 
     labels = butternut_24pop_names, cex = 0.8, pos = 1)

abline(hexp_dist_lm, col = "dodgerblue4")
abline(hexp_wo_wi_lm, col = "darkorchid")

##legends
legend('topleft', legend = hexp_dist_rp, pch = 4, col = "dodgerblue4",
       title = "Regression with WI",
       bty = 'n', border = "black", pt.cex = 1, cex = 0.8)

legend('bottomleft', legend = hexp_wo_wi_dist_rp, pch = 4, 
       col = "darkorchid", title = "Regression without WI",
       bty = 'n', border = "black", pt.cex = 1, cex = 0.8)

legend('bottom', legend = c("New Brunswick", "Ontario", "United States"), pch = 17, col = c("firebrick1", "firebrick4","dodgerblue"))

dev.off()

##create data frame of pvalues r2 for linear and quad models
dist_edge_hexp_df <- matrix(nrow = 2, ncol = 2)
rownames(dist_edge_hexp_df) <- c("linear_hexp_withWI","linear_hexp_withoutWI")
colnames(dist_edge_hexp_df) <- c("Pvalue","R2")
##write values
dist_edge_hexp_df[1,] <- c(hexp_dist_pvalue, hexp_dist_r2)
dist_edge_hexp_df[2,] <- c(hexp_wo_wi_dist_pvalue, hexp_wo_wi_dist_r2)

##write out
write.csv(dist_edge_hexp_df, paste0(dist_edge_path, "hexp_dist_pvalue_r2.csv"))

