##########################
######## Libraries #######
##########################

library(adegenet)
library(poppr)
library(hierfstat)

#####################################
############ Directories ############
#####################################
butternut_drive <- "C:\\Users\\eschumacher\\Documents\\GitHub\\butternut"

##############################################################
#################### Load Genetic Files  #####################
##############################################################
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

##generate poppr file for population size and statistics
butternut_24pop_poppr <- poppr(butternutgen_reorg)

##dist edge df 
dist_edge_df <- read.csv("data_files\\geographic_files\\butternut_dist_edge_df.csv")

############################################################################
############ Calculate Allelic Richness and Create DF  #####################
############################################################################
##calculate allelic richness for each population 
butternut_reorg_allrich <- colSums(allelic.richness(butternutgen_reorg)$Ar)/length(butternutgen_reorg@loc.n.all)

##create df with allelic richness and distance
all_rich_dist_df <- data.frame(dist_edge_df$Dist_To_Edge, butternut_reorg_allrich, dist_edge_df$Col)
colnames(all_rich_dist_df) <- c("Distance", "Allelic_Richness", "Col")
rownames(all_rich_dist_df) <- butternut_24pop_names

##create data frame without Wisconsin populations 
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
pdf("genetic_analyses_results\\linear_allrich_dist_edge.pdf", width = 8, height = 6)

plot(all_rich_dist_df[,2]~all_rich_dist_df[,1], 
     col = as.character(all_rich_dist_df[,3]), pch = 17,
     ylab = "Allelic Richness", xlab = "Distance to Range Edge (km)", 
     cex = (butternut_24pop_poppr[1:24,2]/50), ylim = c(5,10), xlim = c(0,600),
     main = "Allelic Richness Compared with Distance to Range Edge (km)")

text(all_rich_dist_df[,2]~all_rich_dist_df[,1], 
     labels = rownames(all_rich_dist_df), cex = 0.8, pos = 1)

abline(allrich_lm, col = "dodgerblue4", lwd = 1)
abline(allrich_red_lm, col = "darkorchid", lwd = 1)

##legends
legend('topleft', legend = linear_allrich_dist_rp, pch = 17, col = "dodgerblue4",
       title = "Regression with WI",
       bty = 'n', border = "black", pt.cex = 1, cex = 0.8)

legend('bottomleft', legend = linear_allrich_dist_red_rp, pch = 17, 
       col = "darkorchid", title = "Regression without WI",
       bty = 'n', border = "black", pt.cex = 1, cex = 0.8)

legend('bottom', legend = c("New Brunswick", "Ontario", "Quebec", "United States"), pch = 17, 
       col = c("firebrick1", "firebrick4", "lightsalmon", "dodgerblue"))

dev.off()

##create data frame of pvalues r2 for linear and quad models
dist_edge_allrich_df <- matrix(nrow = 2, ncol = 2)
rownames(dist_edge_allrich_df) <- c("linear_allrich_dist_edge",
                                    "linear_allrich_dist_edge_woWI")
colnames(dist_edge_allrich_df) <- c("Pvalue","R2")
##write values
dist_edge_allrich_df[1,] <- c(linear_allrich_dist_pvalue, linear_allrich_dist_r2)
dist_edge_allrich_df[2,] <- c(linear_allrich_dist_red_pvalue, linear_allrich_dist_red_r2)

##write out data frame with r2 and p values
write.csv(dist_edge_allrich_df, "genetic_analyses_results\\pvalue_r2_dist_edge_allrich.csv")

##write out data frame with distance to edge and allelic richness 
write.csv(all_rich_dist_df, "genetic_analyses_results\\allrich_dist_df.csv")
