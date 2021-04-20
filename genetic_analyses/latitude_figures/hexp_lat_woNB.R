##########################
######## Libraries #######
##########################

library(adegenet)
library(poppr)

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

##load relatedness
reorg_relatedness <- read.csv("data_files\\after_reorg\\reorg_relatedness.csv")

##rename genind 
rownames(butternutgen_reorg@tab) <- reorg_relatedness$Ind

##load in mean longitude and latitude by population 
hexp_lat_df <- read.csv("genetic_analyses_results\\hexp_lat_df.csv")

##get population names
butternut_24pop_names <- unique(reorg_relatedness$Pop)

##name populations in genind file 
levels(butternutgen_reorg@pop) <- butternut_24pop_names

##create butternut poppr
butternut_poppr <- poppr(butternutgen_reorg)

##read in maximum and minimum latitudes 
max_min_lonlat_df <- read.csv("data_files\\geographic_files\\max_min_lonlat_df.csv")

##########################################################
#################### Separate NB  ########################
##########################################################
hexp_lat_df_wo_NB <- hexp_lat_df[-c(1:6),]

hexp_lat_df_NB <- hexp_lat_df[c(1:6),]

#############################################
############### Hexp Models  ################
#############################################
####Hexp linear model without NB 
hexp_lat_woNB_lm <- lm(hexp_lat_df[,4]~hexp_lat_df[,2])
hexp_lat_woNB_lm_sum <- summary(hexp_lat_woNB_lm)

##create rp document 
hexp_woNB_r2 <- hexp_lat_woNB_lm_sum$adj.r.squared
hexp_woNB_coef <- hexp_lat_woNB_lm_sum$coefficients
hexp_woNB_pvalue <- hexp_woNB_coef[2,4]
hexp_woNB_rp <- vector('expression',2)
hexp_woNB_rp[1] <- substitute(expression(italic(R)^2 == MYVALUE), 
                         list(MYVALUE = format(hexp_woNB_r2,dig=3)))[2]
hexp_woNB_rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                        list(MYOTHERVALUE = format(hexp_woNB_pvalue, digits = 2)))[2]

###df reduced without wisconsin
hexp_lat_woNB_red_df <- hexp_lat_df_wo_NB[-c(10,16:17),]

####Hexp reduced 
hexp_lat_woNB_red_lm <- lm(hexp_lat_woNB_red_df[,4]~hexp_lat_woNB_red_df[,2])
hexp_lat_woNB_red_lm_sum <- summary(hexp_lat_woNB_red_lm)

##create rp document 
hexp_red_woNB_r2 <- hexp_lat_woNB_red_lm_sum$adj.r.squared
hexp_red_woNB_coef <- hexp_lat_woNB_red_lm_sum$coefficients
hexp_red_woNB_pvalue <- hexp_red_woNB_coef[2,4]
hexp_red_woNB_rp <- vector('expression',2)
hexp_red_woNB_rp[1] <- substitute(expression(italic(R)^2 == MYVALUE), 
                             list(MYVALUE = format(hexp_red_woNB_r2,dig=3)))[2]
hexp_red_woNB_rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                            list(MYOTHERVALUE = format(hexp_red_woNB_pvalue, digits = 2)))[2]

###Linear Regression
pdf("genetic_analyses_results\\hexp_lat_linear_woNB.pdf", width = 8, height = 6)

##plotting code
plot(hexp_lat_df_wo_NB[,4]~hexp_lat_df_wo_NB[,2], col = hexp_lat_df_wo_NB[,3], pch = 17, ylim = c(0.7, 0.9), 
     main = "Expected Heterozygosity Compared to Mean Latitude", ylab = "Expected Heterozygosity",
     xlab = "Mean Latitude", cex = (butternut_poppr[1:24,2]/50))
text(hexp_lat_df_wo_NB[,4]~hexp_lat_df_wo_NB[,2],labels = hexp_lat_df_wo_NB[,1], cex = 0.8, pos = 1)
abline(hexp_lat_woNB_lm, col = "dodgerblue4")
abline(hexp_lat_woNB_red_lm, col = "darkorchid4")
legend('bottom', legend = c("Ontario", "Quebec","United States"), pch = 17, 
       col = c("firebrick4", "lightsalmon", "dodgerblue"))
legend('topleft', legend = hexp_woNB_rp, bty = 'n', border = "black", pt.cex = 1, 
       cex = 1, pch = 17, col = "dodgerblue4",
       title = "With WI populations")
legend('bottomleft', legend = hexp_red_woNB_rp, bty = 'n', border = "black", pt.cex = 1, 
       cex = 1, pch = 17, col = "darkorchid4", title = "Without WI populations")

dev.off()

###############################
##### Quadratic Regression ####
###############################
points_woNB_hexp <- hexp_lat_df_wo_NB[,2]
points2_woNB_hexp <- points_woNB_hexp^2

###Quadratic linear model
quad_hexp_woNB_lm <-lm(hexp_lat_df_wo_NB[,4] ~ 
                                     points_woNB_hexp + points2_woNB_hexp)
##Set up model
quad_hexp_woNB_lm_sum <- summary(quad_hexp_woNB_lm)

#quadratic models
points_hexp_woNB_values <- seq(max_min_lonlat_df[1,3], max_min_lonlat_df[2,3], 0.5)
points_hexp_woNB_counts <- predict(quad_hexp_woNB_lm, list(points_woNB_hexp=points_hexp_woNB_values, 
                                                           points2_woNB_hexp=points_hexp_woNB_values^2))

###points without wisconsin populations
points_hexp_woNB_21 <- hexp_lat_woNB_red_df[,2]
points2_hexp_woNB_21 <- points_hexp_woNB_21^2

######Quadratic 
quad_21pops_woNB_lm <-lm(hexp_lat_woNB_red_df[,4] ~ 
                                       points_hexp_woNB_21 + points2_hexp_woNB_21)
quad_21pops_woNB_lm_sum <- summary(quad_21pops_woNB_lm)

####Quadratic regression
points_hexp_values_woNB_21 <- seq(max_min_lonlat_df[1,3], max_min_lonlat_df[2,3], 0.5)
points_hexp_counts_woNB_21 <- predict(quad_21pops_woNB_lm,list(points_hexp_woNB_21=points_hexp_values_woNB_21, 
                                                                          points2_hexp_woNB_21=points_hexp_values_woNB_21^2))

##now create rp document 
hexp_quad_woNB_r2 <- quad_hexp_woNB_lm_sum$adj.r.squared
hexp_quad_woNB_coef <- quad_hexp_woNB_lm_sum$coefficients
hexp_quad_woNB_pvalue <- hexp_quad_woNB_coef[2,4]
hexp_quad_woNB_rp <- vector('expression',2)
hexp_quad_woNB_rp[1] <- substitute(expression(italic(R)^2 == MYVALUE), 
                              list(MYVALUE = format(hexp_quad_woNB_r2,dig=3)))[2]
hexp_quad_woNB_rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                             list(MYOTHERVALUE = format(hexp_quad_woNB_pvalue, digits = 2)))[2]

##rp document with reduced data frame
hexp_quad_red_woNB_r2 <- quad_21pops_woNB_lm_sum$adj.r.squared
hexp_quad_red_woNB_coef <- quad_21pops_woNB_lm_sum$coefficients
hexp_quad_red_woNB_pvalue <- hexp_quad_red_woNB_coef[2,4]
hexp_quad_red_woNB_rp <- vector('expression',2)
hexp_quad_red_woNB_rp[1] <- substitute(expression(italic(R)^2 == MYVALUE), 
                                  list(MYVALUE = format(hexp_quad_red_woNB_r2,dig=3)))[2]
hexp_quad_red_woNB_rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                                 list(MYOTHERVALUE = format(hexp_quad_red_woNB_pvalue, digits = 2)))[2]

#####plot 
pdf("genetic_analyses_results\\lat_hexp_quad_woNB.pdf", width = 8, height = 6)

plot(hexp_lat_df_wo_NB[,4]~hexp_lat_df_wo_NB[,2], col = hexp_lat_df_wo_NB[,3], pch = 17, 
     main = "Expected Heterozygosity Compared to Mean Latitude with NB", ylab = "Expected Heterozygosity", 
     xlab = "Mean Latitude", cex = (butternut_poppr[1:24,2]/50), ylim = c(0.7, 0.9))
text(hexp_lat_df_wo_NB[,4]~hexp_lat_df_wo_NB[,2], labels = hexp_lat_df_wo_NB[,1], cex = 0.8, pos = 1)
lines(points_hexp_woNB_values, points_hexp_woNB_counts, col = "darkslategray3", lwd = 3)
lines(points_hexp_values_woNB_21, points_hexp_counts_woNB_21, col = "darkseagreen4", lwd = 3)
legend('bottom', legend = c("Ontario", "Quebec","United States"), pch = 17, 
       col = c("firebrick4", "lightsalmon", "dodgerblue"))
legend('topleft', legend = hexp_quad_woNB_rp, bty = 'n', border = "black", pt.cex = 1, cex = 1, pch = 17, 
       col = "darkslategray3", title = "With WI populations")
legend('bottomleft', legend = hexp_quad_red_woNB_rp, bty = 'n', border = "black", pt.cex = 1, cex = 1, pch = 17, 
       col = "darkseagreen4", title = "Without WI populations")

dev.off()

#####write out data frames 
##paste into a df 
hexp_woNB_r2_pvalue_lat <- matrix(nrow = 4, ncol = 2)
rownames(hexp_woNB_r2_pvalue_lat) <- c("Linear_hexp_woNB", "Linear_hexp_woWI_woNB","Quad_Hexp_woNB",
                                       "Quadatric_hexp_woNB_woWI")
colnames(hexp_woNB_r2_pvalue_lat) <- c("Pvalue", "R2")

##Input r2 and pvalues 
hexp_woNB_r2_pvalue_lat[1,] <- c(hexp_woNB_pvalue, hexp_woNB_r2)
hexp_woNB_r2_pvalue_lat[2,] <- c(hexp_red_woNB_pvalue, hexp_red_woNB_r2)
hexp_woNB_r2_pvalue_lat[3,] <- c(hexp_quad_woNB_pvalue, hexp_quad_woNB_r2)
hexp_woNB_r2_pvalue_lat[4,] <- c(hexp_quad_red_woNB_pvalue, hexp_quad_red_woNB_r2)

##write out table
write.csv(hexp_woNB_r2_pvalue_lat,"genetic_analyses_results\\hexp_woNB_r2_pvalue_lat.csv")
