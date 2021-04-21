##########################
######## Libraries #######
##########################

library(adegenet)
library(hierfstat)
library(poppr)

#####################################
############ Directories ############
#####################################
butternut_drive <- "C:\\Users\\eschumacher\\Documents\\GitHub\\butternut"

##############################################################
######################### Load Files  ########################
##############################################################
setwd(butternut_drive)

##load current working reorganized genepop file 
butternutgen_reorg <- read.genepop("data_files\\after_reorg\\reorg_gen_24pop.gen", ncode = 3)

##load relatedness file to name individuals
reorg_relatedness <- read.csv("data_files\\after_reorg\\reorg_relatedness.csv")

###name individuals within genind file
rownames(butternutgen_reorg@tab) <- reorg_relatedness$Ind

##create population name file 
butternut_24pop_names <- unique(reorg_relatedness$Pop)

##name populations within genind file
levels(butternutgen_reorg@pop) <- butternut_24pop_names

##generate poppr file for population size and statistics
butternut_24pop_poppr <- poppr(butternutgen_reorg)

##dist edge df 
dist_edge_df <- read.csv("data_files\\geographic_files\\butternut_dist_edge_df.csv")

######################################################################
##################### Create hexp data frame  ########################
######################################################################
##create hexp df by combining distance to edge, hexp, and color for each pop
hexp_dist_edge_df <- data.frame(dist_edge_df$Dist_To_Edge, as.numeric(butternut_24pop_poppr[1:24,10]), dist_edge_df$Col)
##name rows and columns
colnames(hexp_dist_edge_df) <- c("Distance", "HExp", "Col")
rownames(hexp_dist_edge_df) <- butternut_24pop_names

##remove wisconsin populations 
hexp_dist_edge_red <- hexp_dist_edge_df[-c(16,22,23),]

######################################################################
########################## Linear Models  ############################
######################################################################

###now calculate regressions - linear
hexp_dist_lm <- lm(hexp_dist_edge_df[,2]~hexp_dist_edge_df[,1])
hexp_dist_lm_sum <- summary(hexp_dist_lm)

##linear regression without 
hexp_dist_edge_red_lm <- lm(hexp_dist_edge_red[,2]~hexp_dist_edge_red[,1])
hexp_dist_edge_red_lm_sum <- summary(hexp_dist_edge_red_lm)

##create r2 and p for full df 
hexp_dist_r2 <- hexp_dist_lm_sum$adj.r.squared
hexp_dist_coef <- hexp_dist_lm_sum$coefficients
hexp_dist_pvalue <- hexp_dist_coef[2,4]
hexp_dist_rp <- vector('expression',2)
hexp_dist_rp[1] <- substitute(expression(italic(R)^2 == MYVALUE), 
                                        list(MYVALUE = format(hexp_dist_r2,dig=3)))[2]
hexp_dist_rp[2] <- substitute(expression(italic(p) == MYOTHERVALUE), 
                                        list(MYOTHERVALUE = format(hexp_dist_pvalue, digits = 2)))[2]

##create r2 and p for reduced df - without Wisconsin populations 
hexp_red_dist_r2 <- hexp_dist_edge_red_lm_sum$adj.r.squared
hexp_red_dist_coef <- hexp_dist_edge_red_lm_sum$coefficients
hexp_red_dist_pvalue <- hexp_red_dist_coef[2,4]
hexp_red_dist_rp <- vector('expression',2)
hexp_red_dist_rp[1] <- substitute(expression(italic(R)^2 == MYVALUE), 
                                            list(MYVALUE = format(hexp_red_dist_r2,dig=3)))[2]
hexp_red_dist_rp[2] <- substitute(expression(italic(p) == MYOTHERVALUE), 
                                            list(MYOTHERVALUE = format(hexp_red_dist_pvalue, digits = 2)))[2]

##################################################################
######################## Plotting Models  ########################
##################################################################
##HExp linear model
pdf("genetic_analyses_results\\hexp_dist.pdf", width = 8, height = 6)

plot(hexp_dist_edge_df[,2]~hexp_dist_edge_df[,1], 
     col = as.character(hexp_dist_edge_df[,3]), pch = 17, ylim = c(0.74,0.86),
     ylab = "Expected Heterozygosity", xlab = "Distance to Range Edge (km)", 
     cex = (butternut_24pop_poppr[1:24,2]/50), xlim = c(0,600),
     main = "Expected Heterozygosity Compared with Distance to Range Edge (km)")

text(hexp_dist_edge_df[,2]~hexp_dist_edge_df[,1], 
     labels = butternut_24pop_names, cex = 0.8, pos = 1)

abline(hexp_dist_lm, col = "dodgerblue4")
abline(hexp_dist_edge_red_lm, col = "darkorchid")

##legends
legend('topleft', legend = hexp_dist_rp, pch = 17, col = "dodgerblue4",
       title = "Regression with WI",
       bty = 'n', border = "black", pt.cex = 1, cex = 0.8)

legend('bottomleft', legend = hexp_red_dist_rp, pch = 17, 
       col = "darkorchid", title = "Regression without WI",
       bty = 'n', border = "black", pt.cex = 1, cex = 0.8)

legend('bottomright', legend = c("New Brunswick", "Ontario", "Quebec", "United States"), 
       pch = 17, col = c("firebrick1", "firebrick4", "lightsalmon","dodgerblue"))

dev.off()

##create data frame of pvalues r2 for linear and quad models
dist_edge_hexp_df <- matrix(nrow = 2, ncol = 2)
rownames(dist_edge_hexp_df) <- c("linear_hexp_withWI","linear_hexp_withoutWI")
colnames(dist_edge_hexp_df) <- c("Pvalue","R2")
##write values
dist_edge_hexp_df[1,] <- c(hexp_dist_pvalue, hexp_dist_r2)
dist_edge_hexp_df[2,] <- c(hexp_red_dist_pvalue, hexp_red_dist_r2)

##write out pvalues and r2 relatedness to the relationship
write.csv(dist_edge_hexp_df, "genetic_analyses_results\\hexp_dist_pvalue_r2.csv")

##write out data frame with distance to edge and hexp 
write.csv(hexp_dist_edge_df, "genetic_analyses_results\\hexp_dist_edge_df.csv")
