##########################
######## Libraries #######
##########################

library(adegenet)
library(poppr)

#####################################
############ Directories ############
#####################################
butternut_drive <- "C:\\Users\\eschumacher\\Documents\\GitHub\\butternut"

######################################################
#################### Load Files  #####################
######################################################
##set working directory
setwd(butternut_drive)

##load current working reorganized genepop file 
butternutgen_reorg <- read.genepop("data_files\\after_reorg\\reorg_gen_24pop.gen", ncode = 3)

##load relatedness file to name individuals
reorg_relatedness <- read.csv("data_files\\after_reorg\\reorg_relatedness.csv")

###Name the reorg file
rownames(butternutgen_reorg@tab) <- reorg_relatedness$Ind

##create population name file 
butternut_24pop_names <- unique(reorg_relatedness$Pop)

##Name levels
levels(butternutgen_reorg@pop) <- butternut_24pop_names

##generate a poppr document 
butternut_24pop_poppr <- poppr(butternutgen_reorg)

##load in allelic richness and distance edge document
allrich_dist_edge_df <- read.csv("genetic_analyses_results\\allrich_dist_edge_df.csv")

##############################################################################
######## Run Analyses and Create DF of Just Ontario Quebec and US  ###########
##############################################################################
##remove New Brunswick and create DFs
all_rich_dist_df_wo_NB <- allrich_dist_edge_df[-c(1:6),]

##without wisconsin populations
all_rich_dist_df_wo_NB_red <- all_rich_dist_df_wo_NB[-c(10,16:17),]

######################################################################
########################## Linear Models  ############################
######################################################################
###now calculate regressions - linear
allrich_wo_NB_lm <- lm(all_rich_dist_df_wo_NB[,3]~all_rich_dist_df_wo_NB[,2])
allrich_wo_NB_lm_sum <- summary(allrich_wo_NB_lm)

##Reduced data frame
allrich_red_lm <- lm(all_rich_dist_df_wo_NB_red[,3]~all_rich_dist_df_wo_NB_red[,2])
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
pdf("genetic_analyses_results\\linear_allrich_dist_woNB.pdf", width = 8, height = 6)

plot(all_rich_dist_df_wo_NB[,3]~all_rich_dist_df_wo_NB[,2], 
     col = as.character(all_rich_dist_df_wo_NB[,4]), pch = 17,
     ylab = "Allelic Richness", xlab = "Distance to Range Edge (km)", 
     cex = (butternut_24pop_poppr[7:24,2]/50), ylim = c(5,10), xlim = c(0,600),
     main = "Allelic Richness Compared with Distance to Range Edge (km)")

text(all_rich_dist_df_wo_NB[,3]~all_rich_dist_df_wo_NB[,2], 
     labels = all_rich_dist_df_wo_NB[,1], cex = 0.8, pos = 1)

abline(allrich_wo_NB_lm, col = "dodgerblue4")
abline(allrich_red_lm, col = "darkorchid")

##legends
legend('topleft', legend = linear_allrich_dist_wo_NB_rp, pch = 17, col = "dodgerblue4",
       title = "Regression with WI wo NB",
       bty = 'n', border = "black", pt.cex = 1, cex = 0.8)

legend('bottomleft', legend = linear_allrich_dist_wo_NB_red_rp, pch = 17, 
     col = "darkorchid", title = "Regression without WI wo NB",
     bty = 'n', border = "black", pt.cex = 1, cex = 0.8)

legend('bottom', legend = c("Ontario", "Quebec", "United States"), 
       pch = 17, col = c("firebrick4", "lightsalmon","dodgerblue"))

dev.off()

###write out data frame with p values and r2
##create data frame
dist_edge_allrich_woNB_df <- matrix(nrow = 2, ncol = 2)

##insert values
dist_edge_allrich_woNB_df[1,] <- c(linear_allrich_dist_wo_NB_pvalue,linear_allrich_dist_wo_NB_r2)
dist_edge_allrich_woNB_df[2,] <- c(linear_allrich_dist_wo_NB_red_pvalue, linear_allrich_dist_wo_NB_red_r2)

##name columns and row names 
rownames(dist_edge_allrich_woNB_df) <- c("Distance to Edge Allelic Richness without NB",
                                 "Distance to Edge Allelic Richness without NB and without WI")
colnames(dist_edge_allrich_woNB_df) <- c("P-value","R2")

##write out csv
write.csv(dist_edge_allrich_woNB_df, "genetic_analyses_results\\dist_edge_allrich_woNB_df.csv")


