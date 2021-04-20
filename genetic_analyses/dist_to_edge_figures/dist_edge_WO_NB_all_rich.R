##########################
######## Libraries #######
##########################

library(diveRsity)
library(adegenet)
library(tidyr)
library(hierfstat)
library(poppr)
library(Demerelate)

#####################################
############ Directories ############
#####################################
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
pdf("Graphical_Stat_Results\\PostIndRemoval\\24pop\\Reorg_Results\\GeneticDiversity\\linear_allrich_dist_woNB.pdf", width = 8, height = 6)

plot(all_rich_dist_df[,2]~all_rich_dist_df[,1], 
     col = as.character(all_rich_dist_df[,3]), pch = 17,
     ylab = "Number of Alleles", xlab = "Distance to Range Edge (km)", 
     cex = (butternut_24pop_poppr[7:24,2]/80), ylim = c(5,10), xlim = c(0,600),
     main = "Number of Alleles Compared with Distance to Range Edge (km)")

text(all_rich_dist_df[,2]~all_rich_dist_df[,1], 
     labels = rownames(all_rich_dist_df), cex = 0.8, pos = 1)

abline(allrich_wo_NB_lm, col = "dodgerblue4")
abline(allrich_red_lm, col = "darkorchid")

##legends
legend('topleft', legend = linear_allrich_dist_wo_NB_rp, pch = 17, col = "dodgerblue4",
       title = "Regression with WI wo NB",
       bty = 'n', border = "black", pt.cex = 1, cex = 0.8)

legend('bottomleft', legend = linear_allrich_dist_wo_NB_red_rp, pch = 17, 
     col = "darkorchid", title = "Regression without WI wo NB",
     bty = 'n', border = "black", pt.cex = 1, cex = 0.8)

legend('bottom', legend = c("New Brunswick","Ontario", "United States"), pch = 17, col = c("firebrick1","firebrick4","dodgerblue"))

dev.off()

###########################################
########## Write out R2 and Plot ##########
###########################################
dist_edge_woNB_df <- matrix(nrow = 2, ncol = 2)
dist_edge_woNB_df[1,] <- c(linear_allrich_dist_wo_NB_pvalue,linear_allrich_dist_wo_NB_r2)
dist_edge_woNB_df[2,] <- c(linear_allrich_dist_wo_NB_red_pvalue, linear_allrich_dist_wo_NB_red_r2)

##name columns and row names 
rownames(dist_edge_woNB_df) <- c("Distance to Edge Allelic Richness without NB",
                                 "Distance to Edge Allelic Richness without NB and without WI")
colnames(dist_edge_woNB_df) <- c("P-value","R2")

##write out csv
write.csv(dist_edge_woNB_df, "Graphical_Stat_Results\\PostIndRemoval\\24pop\\Reorg_Results\\GeneticDiversity\\dist_edge_woNB_df.csv")


