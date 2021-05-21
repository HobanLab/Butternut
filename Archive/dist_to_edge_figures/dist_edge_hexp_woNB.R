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
hexp_dist_edge_df <- read.csv("genetic_analyses_results\\hexp_dist_edge_df.csv")

##############################################################################
######## Run Analyses and Create DF of Just Ontario Quebec and US  ###########
##############################################################################
##remove New Brunswick and create DFs
hexp_dist_edge_df_woNB <- hexp_dist_edge_df[-c(1:6),]

##without wisconsin populations
hexp_dist_edge_df_woNB_red <- hexp_dist_edge_df_woNB[-c(10,16:17),]

######################################################################
########################## Linear Models  ############################
######################################################################
###now calculate regressions - linear
hexp_woNB_lm <- lm(hexp_dist_edge_df_woNB[,3]~hexp_dist_edge_df_woNB[,2])
hexp_woNB_lm_sum <- summary(hexp_woNB_lm)

##Reduced data frame
hexp_woNB_red_lm <- lm(hexp_dist_edge_df_woNB_red[,3]~hexp_dist_edge_df_woNB_red[,2])
hexp_woNB_red_lm_sum <- summary(hexp_woNB_red_lm)

##create r2 and p for full df 
hexp_dist_edge_woNB_r2 <- hexp_woNB_lm_sum$adj.r.squared
hexp_dist_edge_woNB_coef <- hexp_woNB_lm_sum$coefficients
hexp_dist_edge_woNB_pvalue <- hexp_dist_edge_woNB_coef[2,4]
hexp_dist_edge_woNB_rp <- vector('expression',2)
hexp_dist_edge_woNB_rp[1] <- substitute(expression(italic(R)^2 == MYVALUE), 
                                              list(MYVALUE = format(hexp_dist_edge_woNB_r2,dig=3)))[2]
hexp_dist_edge_woNB_rp[2] <- substitute(expression(italic(p) == MYOTHERVALUE), 
                                              list(MYOTHERVALUE = format(hexp_dist_edge_woNB_pvalue, digits = 2)))[2]

##create r2 and p for reduced df 
hexp_dist_edge_woNB_red_r2 <- hexp_woNB_red_lm_sum$adj.r.squared
hexp_dist_edge_woNB_red_coef <- hexp_woNB_red_lm_sum$coefficients
hexp_dist_edge_woNB_red_pvalue <- hexp_dist_edge_woNB_red_coef[2,4]
hexp_dist_edge_woNB_red_rp <- vector('expression',2)
hexp_dist_edge_woNB_red_rp[1] <- substitute(expression(italic(R)^2 == MYVALUE), 
                                                  list(MYVALUE = format(hexp_dist_edge_woNB_red_r2,dig=3)))[2]
hexp_dist_edge_woNB_red_rp[2] <- substitute(expression(italic(p) == MYOTHERVALUE), 
                                                  list(MYOTHERVALUE = format(hexp_dist_edge_woNB_red_pvalue, digits = 2)))[2]

##################################################################
################# Plotting Linear Models  ########################
##################################################################

##Allelic Richness Distance
pdf("genetic_analyses_results\\linear_hexp_dist_woNB.pdf", width = 8, height = 6)

plot(hexp_dist_edge_df_woNB[,3]~hexp_dist_edge_df_woNB[,2], 
     col = as.character(hexp_dist_edge_df_woNB[,4]), pch = 17,
     ylab = "Expected Heterozygosity", xlab = "Distance to Range Edge (km)", 
     cex = (butternut_24pop_poppr[7:24,2]/50), ylim = c(0.74,0.89), xlim = c(0,600),
     main = "Expected Heterozygosity Compared with Distance to Range Edge (km)")

text(hexp_dist_edge_df_woNB[,3]~hexp_dist_edge_df_woNB[,2], 
     labels = hexp_dist_edge_df_woNB[,1], cex = 0.8, pos = 1)

abline(hexp_woNB_lm, col = "dodgerblue4")
abline(hexp_woNB_red_lm, col = "darkorchid")

##legends
legend('topleft', legend = hexp_dist_edge_woNB_rp, pch = 17, col = "dodgerblue4",
       title = "Regression with WI wo NB",
       bty = 'n', border = "black", pt.cex = 1, cex = 0.8)

legend('bottomleft', legend = hexp_dist_edge_woNB_red_rp, pch = 17, 
       col = "darkorchid", title = "Regression without WI wo NB",
       bty = 'n', border = "black", pt.cex = 1, cex = 0.8)

legend('bottomright', legend = c("Ontario", "Quebec", "United States"), 
       pch = 17, col = c("firebrick4", "lightsalmon","dodgerblue"))

dev.off()

###write out data frame with p values and r2
##create data frame
dist_edge_hexp_woNB_df <- matrix(nrow = 2, ncol = 2)

##insert values
dist_edge_hexp_woNB_df[1,] <- c(hexp_dist_edge_woNB_pvalue,hexp_dist_edge_woNB_r2)
dist_edge_hexp_woNB_df[2,] <- c(hexp_dist_edge_woNB_red_pvalue, hexp_dist_edge_woNB_red_r2)

##name columns and row names 
rownames(dist_edge_hexp_woNB_df) <- c("Distance to Edge Expected Heterozygosity without NB",
                                         "Distance to Edge Expected Heterozygosity without NB and without WI")
colnames(dist_edge_hexp_woNB_df) <- c("P-value","R2")

##write out csv
write.csv(dist_edge_hexp_woNB_df, "genetic_analyses_results\\dist_edge_hexp_woNB_df.csv")
