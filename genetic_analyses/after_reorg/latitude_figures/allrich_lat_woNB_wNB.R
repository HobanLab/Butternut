##########################
######## Libraries #######
##########################

library(diveRsity)
library(adegenet)
library(tidyr)
library(hierfstat)
library(poppr)
library(PopGenReport)

#####################################
############ Directories ############
#####################################
butternut_drive <- "G:\\My Drive\\Hoban_Lab_Docs\\Projects\\Butternut_JUCI"

#####################################
############# Load Files ############
#####################################
setwd(butternut_drive)

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

##########################################################
#################### Separate NB  ########################
##########################################################
all_rich_lat_df_wo_NB <- all_rich_lat_df[-c(1:6),]

all_rich_lat_df_NB <- all_rich_lat_df[c(1:6),]

##########################################################
########## Linear Modeling without NB  ###################
##########################################################
##Create linear relationship with lat and number of alleles 
allrich_lm <- lm(all_rich_lat_df_wo_NB[,2]~all_rich_lat_df_wo_NB[,1])
allrich_lm_sum <- summary(allrich_lm)

###df reduced without wisconsin
all_rich_lat_df_wo_NB_red <- all_rich_lat_df_wo_NB[-c(10,16:17),]

##Now create linear model 
allrich_red_lm <- lm(all_rich_lat_df_wo_NB_red[,2]~all_rich_lat_df_wo_NB_red[,1])
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
pdf("Graphical_Stat_Results\\PostIndRemoval\\24pop\\Reorg_Results\\GeneticDiversity\\linear_woNB_all_rich_Lat.pdf", width = 8, height = 6)

##first plot
plot(all_rich_lat_df_wo_NB[,2]~all_rich_lat_df_wo_NB[,1], col = all_rich_lat_df_wo_NB[,3], pch = 17, 
     main = "Number of Alleles Compared to Mean Latitude", ylab = "Number of Alleles", 
     xlab = "Mean Latitude", cex = (butternut_poppr[7:24,2]/80), ylim = c(5,10))

##name Ontario and US
text(all_rich_lat_df_wo_NB[,2]~all_rich_lat_df_wo_NB[,1], labels = rownames(all_rich_lat_df_wo_NB), 
     cex = 0.8, pos = 1)

##plot linear models
abline(allrich_lm, col = "dodgerblue4")
abline(allrich_red_lm, col = "darkorchid4")

##create linear model 
legend('topleft', legend = allrich_rp, bty = 'n', border = "black", pt.cex = 1, cex = 0.8, pch = 17, 
       col = "dodgerblue4", title = "With WI")
legend('bottomleft', legend = allrich_red_rp, bty = 'n', border = "black", pt.cex = 1, cex = 0.8, pch = 17, 
       col = "darkorchid4", title = "Without WI")
legend('bottom', legend = c("Ontario", "United States"), pch = 17, col = c("firebrick4","dodgerblue"))
dev.off()

##########################################################
################# Linear Modeling NB  ####################
##########################################################
##Create linear relationship with lat and number of alleles 
allrich_NB_lm <- lm(all_rich_lat_df_NB[,2]~all_rich_lat_df_NB[,1])
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
pdf("Graphical_Stat_Results\\PostIndRemoval\\24pop\\Reorg_Results\\GeneticDiversity\\linear_NB_all_rich_Lat.pdf", width = 8, height = 6)

##first plot
plot(all_rich_lat_df_NB[,2]~all_rich_lat_df_NB[,1], col = all_rich_lat_df_NB[,3], pch = 17, 
     main = "Number of Alleles Compared to Mean Latitude", ylab = "Number of Alleles", 
     xlab = "Mean Latitude", cex = (butternut_poppr[1:6,2]/80), ylim = c(5,10))

##name Ontario and US
text(all_rich_lat_df_NB[,2]~all_rich_lat_df_NB[,1], labels = rownames(all_rich_lat_df_NB), 
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
##Do quadratic
#calculate data points 
wo_NB_points <- all_rich_lat_df_wo_NB[,1]
wo_NB_points2 <- wo_NB_points^2

###Quadratic linear model
wo_NB_quadratic_model_lm <-lm(all_rich_lat_df_wo_NB[,2] ~ 
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
wo_NB_points_red <- all_rich_lat_df_wo_NB_red[,1]
wo_NB_points2_red <- all_rich_lat_df_wo_NB_red[,1]^2

##create quadratic regression documents 
quadratic_wo_NB_red_lm <-lm(all_rich_lat_df_wo_NB_red[,2] ~ 
                              wo_NB_points_red + wo_NB_points2_red)
quadratic_wo_NB_red_lm_sum <- summary(quadratic_wo_NB_red_lm)

##values and pointsance
wo_NB_points_values_red <- seq(min(all_rich_lat_df_wo_NB_red[,1]), 
                               max(all_rich_lat_df_wo_NB_red[,1]), 0.5)
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
pdf("Graphical_Stat_Results\\PostIndRemoval\\24pop\\Reorg_Results\\GeneticDiversity\\wo_NB_allrich_quad.pdf", width = 8, height = 6)

###Plot same thing 
plot(all_rich_lat_df_wo_NB[,2]~all_rich_lat_df_wo_NB[,1], col = all_rich_lat_df_wo_NB[,3], 
     pch = 17, main = "Number of Alleles Compared to Mean Latitude", ylab = "Number of Alleles", 
     xlab = "Mean Latitude", xlim = c(35, 47),
     cex = (butternut_poppr[7:24,2]/100), ylim = c(5,10))
points(all_rich_lat_df_NB[,2]~all_rich_lat_df_NB[,1], col = all_rich_lat_df_NB[,3],
       pch = 17)

##Add text
text(all_rich_lat_df[,2]~all_rich_lat_df[,1], labels = butternut_24pop_names, cex = 0.8, pos = 1)

##Now plot lines 
lines(wo_NB_points_values, wo_NB_points_counts, col = "darkslategray3", lwd = 3)
lines(wo_NB_points_values_red, wo_NB_points_counts_red, col = "darkseagreen4", lwd = 3)

##Create legends 
legend('bottom', legend = c("New Brunswick", "Ontario", "United States"), pch = 17, 
       col = c("firebrick1", "firebrick4","dodgerblue"))
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
write.csv(allrich_r2_pvalue_quad_linear_wo_NB_df,"Graphical_Stat_Results\\PostIndRemoval\\24pop\\Reorg_Results\\GeneticDiversity\\allrich_r2_pvalue_quad_linear_wo_NB_df.csv")



