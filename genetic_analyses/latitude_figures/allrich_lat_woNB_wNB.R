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

#####################################
############# Load Files ############
#####################################
##set working directory for files 
setwd(butternut_drive)

##load current working reorganized genepop file 
butternutgen_reorg <- read.genepop("data_files\\after_reorg\\reorg_gen_24pop.gen", ncode = 3)

##load relatedness
reorg_relatedness <- read.csv("data_files\\after_reorg\\reorg_relatedness.csv")

##rename genind 
rownames(butternutgen_reorg@tab) <- reorg_relatedness$Ind

##load in mean longitude and latitude by population 
all_rich_lat_df <- read.csv("genetic_analyses_results\\all_rich_lat_df.csv")

##get pop names 
butternut_24pop_names <- all_rich_lat_df[,1]

##name populations in genind file 
levels(butternutgen_reorg@pop) <- butternut_24pop_names

##create butternut poppr
butternut_poppr <- poppr(butternutgen_reorg)

##########################################################
#################### Separate NB  ########################
##########################################################
all_rich_lat_df_wo_NB <- all_rich_lat_df[-c(1:6),]

all_rich_lat_df_NB <- all_rich_lat_df[c(1:6),]

##########################################################
########## Linear Modeling without NB  ###################
##########################################################
##Create linear relationship with lat and number of alleles 
allrich_lm <- lm(all_rich_lat_df_wo_NB[,3]~all_rich_lat_df_wo_NB[,2])
allrich_lm_sum <- summary(allrich_lm)

###df reduced without wisconsin
all_rich_lat_df_wo_NB_red <- all_rich_lat_df_wo_NB[-c(10,16:17),]

##Now create linear model 
allrich_red_lm <- lm(all_rich_lat_df_wo_NB_red[,3]~all_rich_lat_df_wo_NB_red[,2])
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
pdf("genetic_analyses_results\\linear_woNB_all_rich_Lat.pdf", width = 8, height = 6)

##first plot
plot(all_rich_lat_df_wo_NB[,3]~all_rich_lat_df_wo_NB[,2], col = all_rich_lat_df_wo_NB[,4], pch = 17, 
     main = "Allelic Richness Compared to Mean Latitude", ylab = "Allelic Richness", 
     xlab = "Mean Latitude", cex = (butternut_poppr[7:24,2]/50), ylim = c(5,10))

##name Ontario and US
text(all_rich_lat_df_wo_NB[,3]~all_rich_lat_df_wo_NB[,2], labels = all_rich_lat_df_wo_NB[,1], 
     cex = 0.8, pos = 1)

##plot linear models
abline(allrich_lm, col = "dodgerblue4")
abline(allrich_red_lm, col = "darkorchid4")

##create linear model 
legend('topleft', legend = allrich_rp, bty = 'n', border = "black", pt.cex = 1, cex = 0.8, pch = 17, 
       col = "dodgerblue4", title = "With WI")
legend('bottomleft', legend = allrich_red_rp, bty = 'n', border = "black", pt.cex = 1, cex = 0.8, pch = 17, 
       col = "darkorchid4", title = "Without WI")
legend('bottom', legend = c("Ontario", "Quebec", "United States"), pch = 17, col = c("firebrick4", "lightsalmon",
                                                                                     "dodgerblue"))
dev.off()

##########################################################
################# Linear Modeling NB  ####################
##########################################################
##Create linear relationship with lat and number of alleles 
allrich_NB_lm <- lm(all_rich_lat_df_NB[,3]~all_rich_lat_df_NB[,2])
allrich_NB_lm_sum <- summary(allrich_NB_lm)

##create r2 and p values
allrich_NB_r2 <- allrich_NB_lm_sum$adj.r.squared
allrich_NB_coef <- allrich_NB_lm_sum$coefficients
allrich_NB_pvalue <- allrich_NB_coef[2,4]
allrich_NB_rp <- vector('expression',2)
allrich_NB_rp[1] <- substitute(expression(italic(R)^2 == MYVALUE), 
                            list(MYVALUE = format(allrich_NB_r2,dig=3)))[2]
allrich_NB_rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                           list(MYOTHERVALUE = format(allrich_NB_pvalue, digits = 2)))[2]

##Allelic Richness Linear
pdf("genetic_analyses_results\\linear_NB_all_rich_Lat.pdf", width = 8, height = 6)

##first plot
plot(all_rich_lat_df_NB[,3]~all_rich_lat_df_NB[,2], col = all_rich_lat_df_NB[,4], pch = 17, 
     main = "Number of Alleles Compared to Mean Latitude", ylab = "Number of Alleles", 
     xlab = "Mean Latitude", cex = (butternut_poppr[1:6,2]/80), ylim = c(5,10))

##name Ontario and US
text(all_rich_lat_df_NB[,3]~all_rich_lat_df_NB[,2], labels = all_rich_lat_df_NB[,1], 
     cex = 0.8, pos = 1)

##plot linear models
abline(allrich_NB_lm, col = "dodgerblue4")

##create linear model 
legend('topleft', legend = allrich_NB_rp, bty = 'n', border = "black", pt.cex = 1, cex = 0.8, pch = 17, 
       col = "dodgerblue4")

legend('bottom', legend = c("New Brunswick"), pch = 17, col = c("firebrick1"))
dev.off()

##########################################################
################### Quadratic wo NB ######################
##########################################################
#calculate data points 
wo_NB_points <- all_rich_lat_df_wo_NB[,2]
wo_NB_points2 <- wo_NB_points^2

###Quadratic linear model
wo_NB_quadratic_model_lm <-lm(all_rich_lat_df_wo_NB[,3] ~ 
                                wo_NB_points + wo_NB_points2)
##Calculate summary documents 
wo_NB_quadratic_model_lm_sum <- summary(wo_NB_quadratic_model_lm)

#calculate model out of points 
wo_NB_points_values <- seq(min(all_rich_lat_df_wo_NB$Mean_Lat), 
                           max(all_rich_lat_df_wo_NB$Mean_Lat), 0.5)
wo_NB_points_counts <- predict(wo_NB_quadratic_model_lm,
                               list(wo_NB_points = wo_NB_points_values, 
                                    wo_NB_points2 = wo_NB_points_values^2))

###########without Wisconsin populations
wo_NB_points_red <- all_rich_lat_df_wo_NB_red[,2]
wo_NB_points2_red <- wo_NB_points_red^2

##create quadratic regression documents 
quadratic_wo_NB_red_lm <-lm(all_rich_lat_df_wo_NB_red[,3] ~ 
                              wo_NB_points_red + wo_NB_points2_red)
quadratic_wo_NB_red_lm_sum <- summary(quadratic_wo_NB_red_lm)

##values and pointsance
wo_NB_points_values_red <- seq(min(all_rich_lat_df_wo_NB_red[,2]), 
                               max(all_rich_lat_df_wo_NB_red[,2]), 0.5)
wo_NB_points_counts_red <- predict(quadratic_wo_NB_red_lm,list(wo_NB_points_red = wo_NB_points_values_red, 
                                                               wo_NB_points2_red = wo_NB_points_values_red^2))

##create r2 and p values
wo_NB_allrich_quad_r2 <- wo_NB_quadratic_model_lm_sum$adj.r.squared
wo_NB_allrich_quad_coef <- wo_NB_quadratic_model_lm_sum$coefficients
wo_NB_allrich_quad_pvalue <- wo_NB_allrich_quad_coef[2,4]
wo_NB_allrich_quad_rp <- vector('expression',2)
wo_NB_allrich_quad_rp[1] <- substitute(expression(italic(R)^2 == MYVALUE), 
                                 list(MYVALUE = format(wo_NB_allrich_quad_r2,dig=3)))[2]
wo_NB_allrich_quad_rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                                list(MYOTHERVALUE = format(wo_NB_allrich_quad_pvalue, digits = 2)))[2]

##create r2 and p values for reduced df 
wo_NB_allrich_quad_red_r2 <- quadratic_wo_NB_red_lm_sum$adj.r.squared
wo_NB_allrich_quad_red_coef <- quadratic_wo_NB_red_lm_sum$coefficients
wo_NB_allrich_quad_red_pvalue <- wo_NB_allrich_quad_red_coef[2,4]
wo_NB_allrich_quad_red_rp <- vector('expression',2)
wo_NB_allrich_quad_red_rp[1] <- substitute(expression(italic(R)^2 == MYVALUE), 
                                     list(MYVALUE = format(wo_NB_allrich_quad_red_r2,dig=3)))[2]
wo_NB_allrich_quad_red_rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                                    list(MYOTHERVALUE = format(wo_NB_allrich_quad_red_pvalue, digits = 2)))[2]

########Now Plot quadatric
pdf("genetic_analyses_results\\wo_NB_allrich_quad.pdf", width = 8, height = 6)

###Plot same thing 
plot(all_rich_lat_df_wo_NB[,3]~all_rich_lat_df_wo_NB[,2], col = all_rich_lat_df_wo_NB[,4], 
     pch = 17, main = "Allelic Richness Compared to Mean Latitude", ylab = "Allelic Richness", 
     xlab = "Mean Latitude", xlim = c(35, 47),
     cex = (butternut_poppr[7:24,2]/50), ylim = c(5,10))

##Add text
text(all_rich_lat_df_wo_NB[,3]~all_rich_lat_df_wo_NB[,2], labels = all_rich_lat_df_wo_NB[,1], cex = 0.8, pos = 1)

##Now plot lines 
lines(wo_NB_points_values, wo_NB_points_counts, col = "darkslategray3", lwd = 3)
lines(wo_NB_points_values_red, wo_NB_points_counts_red, col = "darkseagreen4", lwd = 3)

##Create legends 
legend('bottom', legend = c("Ontario", "Quebec", "United States"), pch = 17, 
       col = c("firebrick4", "lightsalmon","dodgerblue"))
legend('topleft', legend = wo_NB_allrich_quad_rp, bty = 'n', border = "black", pt.cex = 1, cex = 0.8, pch = 17, 
       col = "darkslategray3", title = "Quad Regression wo NB")
legend('bottomleft', legend = wo_NB_allrich_quad_red_rp, bty = 'n', border = "black", pt.cex = 1, cex = 0.8, 
       pch = 17, col = "darkseagreen4",
       title = "Quad Regression wo NB and WI")
dev.off()

##paste into a df 
allrich_r2_pvalue_quad_linear_wo_NB_df <- matrix(nrow = 5, ncol = 2)
rownames(allrich_r2_pvalue_quad_linear_wo_NB_df) <- c("Linear_allrich_woNB", "Quadatric_allrich_woNB",
                                                      "Linear_allrich_woNB_woWI","Quadatric_allrich_woNB_woWI",
                                                      "Linear_allrich_NB")
colnames(allrich_r2_pvalue_quad_linear_wo_NB_df) <- c("Pvalue", "R2")

##Input r2 and pvalues 
allrich_r2_pvalue_quad_linear_wo_NB_df[1,] <- c(allrich_pvalue, allrich_r2)
allrich_r2_pvalue_quad_linear_wo_NB_df[2,] <- c(wo_NB_allrich_quad_pvalue, wo_NB_allrich_quad_r2)
allrich_r2_pvalue_quad_linear_wo_NB_df[3,] <- c(allrich_red_pvalue, allrich_red_r2)
allrich_r2_pvalue_quad_linear_wo_NB_df[4,] <- c(wo_NB_allrich_quad_red_pvalue, wo_NB_allrich_quad_red_r2)
allrich_r2_pvalue_quad_linear_wo_NB_df[5,] <- c(allrich_NB_pvalue,allrich_NB_r2)

##write out table
write.csv(allrich_r2_pvalue_quad_linear_wo_NB_df,"genetic_analyses_results\\allrich_r2_pvalue_quad_linear_wo_NB_df.csv")



