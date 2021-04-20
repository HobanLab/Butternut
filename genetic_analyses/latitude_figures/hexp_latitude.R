##########################
######## Libraries #######
##########################

library(adegenet)
library(poppr)

#####################################
############ Directories ############
#####################################
##working directory to load files in 
butternut_drive <- "C:\\Users\\eschumacher\\Documents\\GitHub\\butternut"

#####################################
############# Load Files ############
#####################################
setwd(butternut_drive)

##load current working reorganized genepop file 
butternutgen_reorg <- read.genepop("data_files\\after_reorg\\reorg_gen_24pop.gen", ncode = 3)

##load relatedness
reorg_relatedness <- read.csv("data_files\\after_reorg\\reorg_relatedness.csv")

##rename genind 
rownames(butternutgen_reorg@tab) <- reorg_relatedness$Ind

##load in mean longitude and latitude by population 
butternut_mean_lat_allrich <- read.csv("genetic_analyses_results\\all_rich_lat_df.csv")

##get pop names 
butternut_24pop_names <- unique(reorg_relatedness$Pop)

##name populations in genind file 
levels(butternutgen_reorg@pop) <- butternut_24pop_names

##create butternut poppr
butternut_poppr <- poppr(butternutgen_reorg)

##read in maximum and minimum latitudes 
max_min_lonlat_df <- read.csv("data_files\\geographic_files\\max_min_lonlat_df.csv")

#############################################
############### Hexp Models  ################
#############################################
##################Hexp graph
hexp_lat_df <- cbind(butternut_mean_lat_allrich[,c(1:2,4)], butternut_poppr[1:24,10])
colnames(hexp_lat_df) <- c("Mean Lat", "HExp")
rownames(hexp_lat_df) <- butternut_24pop_names

####Hexp
hexp_lat_lm <- lm(hexp_lat_df[,4]~hexp_lat_df[,2])
hexp_lat_lm_sum <- summary(hexp_lat_lm)

##create rp document 
hexp_r2 <- hexp_lat_lm_sum$adj.r.squared
hexp_coef <- hexp_lat_lm_sum$coefficients
hexp_pvalue <- hexp_coef[2,4]
hexp_rp <- vector('expression',2)
hexp_rp[1] <- substitute(expression(italic(R)^2 == MYVALUE), 
                            list(MYVALUE = format(hexp_r2,dig=3)))[2]
hexp_rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                           list(MYOTHERVALUE = format(hexp_pvalue, digits = 2)))[2]

###df reduced without wisconsin
hexp_lat_red_df <- hexp_lat_df[-c(16,22,23),]

####Hexp reduced 
hexp_lat_red_lm <- lm(hexp_lat_red_df[,4]~hexp_lat_red_df[,2])
hexp_lat_red_lm_sum <- summary(hexp_lat_red_lm)

##create rp document 
hexp_red_r2 <- hexp_lat_red_lm_sum$adj.r.squared
hexp_red_coef <- hexp_lat_red_lm_sum$coefficients
hexp_red_pvalue <- hexp_red_coef[2,4]
hexp_red_rp <- vector('expression',2)
hexp_red_rp[1] <- substitute(expression(italic(R)^2 == MYVALUE), 
                         list(MYVALUE = format(hexp_red_r2,dig=3)))[2]
hexp_red_rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                        list(MYOTHERVALUE = format(hexp_red_pvalue, digits = 2)))[2]

###Linear Regression
pdf("genetic_analyses_results\\hexp_lat_linear.pdf", width = 8, height = 6)

##plotting code
plot(hexp_lat_df[,4]~hexp_lat_df[,2], col = hexp_lat_df[,3], pch = 17, ylim = c(0.7, 0.9), 
     main = "Expected Heterozygosity Compared to Mean Latitude", ylab = "Expected Heterozygosity",
     xlab = "Mean Latitude", cex = (butternut_poppr[1:24,2]/50))
text(hexp_lat_df[,4]~hexp_lat_df[,2],labels = butternut_24pop_names, cex = 0.8, pos = 1)
abline(hexp_lat_lm, col = "dodgerblue4")
abline(hexp_lat_red_lm, col = "darkorchid4")
legend('bottom', legend = c("New Brunswick", "Ontario", "Quebec","United States"), pch = 17, 
       col = c("firebrick1", "firebrick4", "lightsalmon", "dodgerblue"))
legend('topleft', legend = hexp_rp, bty = 'n', border = "black", pt.cex = 1, cex = 1, pch = 17, col = "dodgerblue4",
       title = "With WI populations")
legend('bottomleft', legend = hexp_red_rp, bty = 'n', border = "black", pt.cex = 1, 
       cex = 1, pch = 17, col = "darkorchid4", title = "Without WI populations")

dev.off()

##Quadratic Regression
####points
points_hexp <- hexp_lat_df[,2]
points2_hexp <- points_hexp^2

###Quadratic linear model
quadratic_model_hexp_lm <-lm(hexp_lat_df[,4] ~ 
                               points_hexp + points2_hexp)
##Set up model
quadratic_model_hexp_lm_sum <- summary(quadratic_model_hexp_lm)

#pointsance models
points_hexp_values <- seq(max_min_lonlat_df[1,3], max_min_lonlat_df[2,3], 0.5)
points_hexp_counts <- predict(quadratic_model_hexp_lm,list(points_hexp=points_hexp_values, 
                                                           points2_hexp=points_hexp_values^2))

###pointsance without wisconsin populations
points_hexp_21 <- hexp_lat_red_df[,2]
points2_hexp_21 <- hexp_lat_red_df[,2]^2

######Quadratic 
quadratic_21pops_model_lm <-lm(hexp_lat_red_df[,4] ~ 
                                 points_hexp_21 + points2_hexp_21)
quadratic_21pops_model_lm_sum <- summary(quadratic_21pops_model_lm)

####Quadratic regression
points_hexp_values_21 <- seq(max_min_lonlat_df[1,3], max_min_lonlat_df[2,3], 0.5)
points_hexp_counts_21 <- predict(quadratic_21pops_model_lm,list(points_hexp_21=points_hexp_values_21, 
                                                                points2_hexp_21=points_hexp_values_21^2))

##now create rp document 
hexp_quad_r2 <- quadratic_model_hexp_lm_sum$adj.r.squared
hexp_quad_coef <- quadratic_model_hexp_lm_sum$coefficients
hexp_quad_pvalue <- hexp_quad_coef[2,4]
hexp_quad_rp <- vector('expression',2)
hexp_quad_rp[1] <- substitute(expression(italic(R)^2 == MYVALUE), 
                                 list(MYVALUE = format(hexp_quad_r2,dig=3)))[2]
hexp_quad_rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                                list(MYOTHERVALUE = format(hexp_quad_pvalue, digits = 2)))[2]

##rp document with reduced data frame
hexp_quad_red_r2 <- quadratic_21pops_model_lm_sum$adj.r.squared
hexp_quad_red_coef <- quadratic_21pops_model_lm_sum$coefficients
hexp_quad_red_pvalue <- hexp_quad_red_coef[2,4]
hexp_quad_red_rp <- vector('expression',2)
hexp_quad_red_rp[1] <- substitute(expression(italic(R)^2 == MYVALUE), 
                              list(MYVALUE = format(hexp_quad_red_r2,dig=3)))[2]
hexp_quad_red_rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                             list(MYOTHERVALUE = format(hexp_quad_red_pvalue, digits = 2)))[2]

#####plot 
pdf("genetic_analyses_results\\lat_hexp_quad.pdf", width = 8, height = 6)

plot(hexp_lat_df[,4]~hexp_lat_df[,2], col = hexp_lat_df[,3], pch = 17, 
     main = "Expected Heterozygosity Compared to Mean Latitude", ylab = "Expected Heterozygosity", 
     xlab = "Mean Latitude", cex = (butternut_poppr[1:24,2]/50), ylim = c(0.7, 0.9))
text(hexp_lat_df[,4]~hexp_lat_df[,2], labels = hexp_lat_df[,1], cex = 0.8, pos = 1)
lines(points_hexp_values, points_hexp_counts, col = "darkslategray3", lwd = 3)
lines(points_hexp_values_21, points_hexp_counts_21, col = "darkseagreen4", lwd = 3)
legend('bottom', legend = c("New Brunswick", "Ontario", "Quebec","United States"), pch = 17, 
       col = c("firebrick1", "firebrick4", "lightsalmon", "dodgerblue"))
legend('topleft', legend = hexp_quad_rp, bty = 'n', border = "black", pt.cex = 1, cex = 1, pch = 17, 
       col = "darkslategray3", title = "With WI populations")
legend('bottomleft', legend = hexp_quad_red_rp, bty = 'n', border = "black", pt.cex = 1, cex = 1, pch = 17, 
       col = "darkseagreen4", title = "Without WI populations")

dev.off()

#####write out data frames 
##paste into a df 
hexp_r2_pvalue_lat <- matrix(nrow = 4, ncol = 2)
rownames(hexp_r2_pvalue_lat) <- c("Linear_hexp", "Linear_hexp_woWI","Quad_Hexp","Quadatric_hexp_woWI")
colnames(hexp_r2_pvalue_lat) <- c("Pvalue", "R2")

##Input r2 and pvalues 
hexp_r2_pvalue_lat[1,] <- c(hexp_pvalue, hexp_r2)
hexp_r2_pvalue_lat[2,] <- c(hexp_red_pvalue, hexp_red_r2)
hexp_r2_pvalue_lat[3,] <- c(hexp_quad_pvalue, hexp_quad_r2)
hexp_r2_pvalue_lat[4,] <- c(hexp_red_pvalue, hexp_quad_red_r2)

##write out p values and r2
write.csv(hexp_r2_pvalue_lat,"genetic_analyses_results\\hexp_r2_pvalue_lat.csv")

##write out heterozygosity 
write.csv(hexp_lat_df, "genetic_analyses_results\\hexp_lat_df.csv")
