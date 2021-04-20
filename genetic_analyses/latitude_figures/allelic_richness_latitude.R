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

#####################################
############# Load Files ############
#####################################
setwd(butternut_drive)

##load current working reorganized genepop file 
butternutgen_reorg <- read.genepop("data_files\\after_reorg\\reorg_gen_24pop.gen", ncode = 3)

##load relatedness document to name individuals 
reorg_relatedness <- read.csv("data_files\\after_reorg\\reorg_relatedness.csv")

##name individuals in the genind file  
rownames(butternutgen_reorg@tab) <- reorg_relatedness$Ind

##load in mean longitude and latitude by population 
butternut_mean_lonlat <- read.csv("data_files\\geographic_files\\butternut_coord_df.csv")

##get pop names 
butternut_24pop_names <- unique(reorg_relatedness$Pop)

##name butternut lon and lat 
butternut_mean_lonlat[,1] <- butternut_24pop_names

##name populations in genind file 
levels(butternutgen_reorg@pop) <- butternut_24pop_names

##create butternut poppr
butternut_poppr <- poppr(butternutgen_reorg)

##load min and maximum latitude and longitude document 
butternut_lonlat_max_min_df <- read.csv("data_files\\geographic_files\\max_min_lonlat_df.csv")

##########################################################
################### Linear Modeling ######################
##########################################################
##Calculate Allelic Richness
reorg_allrich <- colSums(allelic.richness(butternutgen_reorg)$Ar)/length(butternutgen_reorg@loc.n.all)

##create df with allelic richness and latitude
all_rich_lat_df <- data.frame(butternut_mean_lonlat[,3], reorg_allrich)
colnames(all_rich_lat_df) <- c("Mean_Lat", "Allelic_Richness")
rownames(all_rich_lat_df) <- butternut_24pop_names
all_rich_lat_df$Color <- NA
all_rich_lat_df[1:6,3] <- "firebrick1"
all_rich_lat_df[c(8,11),3] <- "lightsalmon"
all_rich_lat_df[c(7,9:10),3] <- "firebrick4"
all_rich_lat_df[12:24,3] <- "dodgerblue"

##Create linear relationship with lat and number of alleles 
allrich_lm <- lm(all_rich_lat_df[,2]~all_rich_lat_df[,1])
allrich_lm_sum <- summary(allrich_lm)

###df reduced without wisconsin
allrich_lat_red_df <- all_rich_lat_df[-c(16,22,23),]

##Now create linear model 
allrich_red_lm <- lm(allrich_lat_red_df[,2]~allrich_lat_red_df[,1])
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
pdf("genetic_analyses_results\\all_rich_lat_linear.pdf", width = 8, height = 6)

##create plot
plot(all_rich_lat_df[,2]~all_rich_lat_df[,1], col = all_rich_lat_df[,3], pch = 17, 
     main = "Number of Alleles Compared to Mean Latitude", ylab = "Number of Alleles", 
     xlab = "Mean Latitude", cex = (butternut_poppr[1:24,2]/50), ylim = c(5,10))
##write text
text(all_rich_lat_df[,2]~all_rich_lat_df[,1], labels = butternut_24pop_names, cex = 0.8, pos = 1)
##write a linear model
abline(allrich_lm, col = "dodgerblue4")
abline(allrich_red_lm, col = "darkorchid4")
##Write legends
legend('topleft', legend = allrich_rp, bty = 'n', border = "black", 
       pt.cex = 1, cex = 0.8, pch = 17, col = "dodgerblue4",
       title = "With WI Populations")
legend('bottomleft', legend = allrich_red_rp, bty = 'n', border = "black", 
       pt.cex = 1, cex = 0.8, pch = 17, col = "darkorchid4",
       title = "Without WI Populations")
legend('bottom', legend = c("New Brunswick", "Ontario", "Quebec","United States"), 
       pch = 17, col = c("firebrick1", "firebrick4","lightsalmon","dodgerblue"))
dev.off()

##########################################################
################# Quadratic Modeling #####################
##########################################################
#calculate data points 
points <- all_rich_lat_df[,1]
points2 <- points^2

###Quadratic linear model
quadratic_model_lm <-lm(all_rich_lat_df[,2] ~ 
                          points + points2)
##Calculate summary documents 
quadratic_model_lm_sum <- summary(quadratic_model_lm)

#calculate model out of points 
points_values <- seq(butternut_lonlat_max_min_df[1,3], butternut_lonlat_max_min_df[2,3], 0.5)
points_counts <- predict(quadratic_model_lm,list(points=points_values, points2=points_values^2))

###pointsance without wisconsin populations
points_21 <- allrich_lat_red_df[,1]
points2_21 <- allrich_lat_red_df[,1]^2

##create quadratic regression documents 
quadratic_21pops_model_lm <-lm(allrich_lat_red_df[,2] ~ 
                                 points_21 + points2_21)
quadratic_21pops_model_lm_sum <- summary(quadratic_21pops_model_lm)

##values and pointsance
points_values_21 <- seq(butternut_lonlat_max_min_df[1,3], butternut_lonlat_max_min_df[2,3], 0.5)
points_counts_21 <- predict(quadratic_21pops_model_lm,list(points_21=points_values_21, points2_21=points_values_21^2))

##create r2 and p values
allrich_quad_r2 <- quadratic_model_lm_sum$adj.r.squared
allrich_quad_coef <- quadratic_model_lm_sum$coefficients
allrich_quad_pvalue <- allrich_quad_coef[2,4]
allrich_quad_rp <- vector('expression',2)
allrich_quad_rp[1] <- substitute(expression(italic(R)^2 == MYVALUE), 
                            list(MYVALUE = format(allrich_quad_r2,dig=3)))[2]
allrich_quad_rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                           list(MYOTHERVALUE = format(allrich_quad_pvalue, digits = 2)))[2]

##create r2 and p values for reduced df 
allrich_quad_red_r2 <- quadratic_21pops_model_lm_sum$adj.r.squared
allrich_quad_red_coef <- quadratic_21pops_model_lm_sum$coefficients
allrich_quad_red_pvalue <- allrich_quad_red_coef[2,4]
allrich_quad_red_rp <- vector('expression',2)
allrich_quad_red_rp[1] <- substitute(expression(italic(R)^2 == MYVALUE), 
                                list(MYVALUE = format(allrich_quad_red_r2,dig=3)))[2]
allrich_quad_red_rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                               list(MYOTHERVALUE = format(allrich_quad_red_pvalue, digits = 2)))[2]

##Now Plot 
pdf("genetic_analyses_results\\all_rich_lat_quad.pdf", width = 8, height = 6)

##start plot 
plot(all_rich_lat_df[,2]~all_rich_lat_df[,1], col = all_rich_lat_df[,3], 
     pch = 17, main = "Allelic Richness Compared to Mean Latitude", ylab = "Allelic Richness", 
     xlab = "Mean Latitude", 
     cex = (butternut_poppr[1:24,2]/50), ylim = c(5,10))
##label text
text(all_rich_lat_df[,2]~all_rich_lat_df[,1], labels = butternut_24pop_names, cex = 0.8, pos = 1)
##draw quadratic lines
lines(points_values, points_counts, col = "darkslategray3", lwd = 3)
lines(points_values_21, points_counts_21, col = "darkseagreen4", lwd = 3)
##draw legends
legend('bottom', legend = c("New Brunswick", "Ontario", "Quebec", "United States"), pch = 17, 
       col = c("firebrick1", "firebrick4", "lightsalmon","dodgerblue"))
legend('topleft', legend = allrich_quad_rp, bty = 'n', border = "black", 
       pt.cex = 1, cex = 0.8, pch = 17, col = "darkslategray3",
       title = "With WI Populations")
legend('bottomleft', legend = allrich_quad_red_rp, bty = 'n', border = "black", 
       pt.cex = 1, cex = 0.8, pch = 17, col = "darkseagreen4",
       title = "Without WI Populations")
dev.off()

##create df of r2 and pvalues 
lat_all_rich_df <- matrix(nrow = 4, ncol = 2)
rownames(lat_all_rich_df) <- c("Linear_Lat_AllRich", "Linear_Lat_AllRich_woWI",
                               "Quadratic_Lat_AllRich", "Quadratic_Lat_AllRich_woWI")
colnames(lat_all_rich_df) <- c("Pvalues","R2")
##add values
lat_all_rich_df[1,] <- c(allrich_pvalue,allrich_r2)
lat_all_rich_df[2,] <- c(allrich_red_pvalue, allrich_red_r2)
lat_all_rich_df[3,] <- c(allrich_quad_pvalue, allrich_quad_r2)
lat_all_rich_df[4,] <- c(allrich_quad_red_pvalue, allrich_quad_red_r2)

##write out csv
write.csv(lat_all_rich_df, "genetic_analyses_results\\r2_pvalue_lat_all_rich_df.csv")
