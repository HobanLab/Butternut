
##########################################################
################# Quadratic Modeling #####################
##########################################################
#We have a nested loop below.  The outer loop goes over the statistics, alleles and heterozygosity
#The inner loop goes over the four subsets of population groupings
#These are the subsets: all pops, all minus WI, all minus NB, all minus (NB+WI)
which_pops<-list(1:24,c(1:15,17:21,24),7:24,c(7:15,17:21,24))
#The following are options used in the loop of 4
line_colors<-c("dodgerblue4","darkorchid4","dodgerblue4","darkorchid4")
loc_legend<-c("topleft","bottomleft","topleft","bottomleft")
title_legend<-c("With WI Populations", "Without WI Populations","With WI Populations", "Without WI Populations")

reorg_allrich <- colSums(allelic.richness(butternutgen_reorg)$Ar)/length(butternutgen_reorg@loc.n.all)
##create df with allelic richness and latitude
all_rich_lat_df <- data.frame(butternut_mean_lonlat[,3], reorg_allrich)
hexp_lat_df <- read.csv("genetic_analyses_results\\hexp_lat_df.csv")[,c(2,4)]

#The loop will go over the two statistics, 1 = allelic richness, 2 = heterozygosity
for (j in 1:2){
  
  #This will set y limits and axis titles to match the statistic
  if (j==1) {stat_lat_df<-all_rich_lat_df; y_low<-5; y_high<-10; stat_name<-"Allelic_Richness"}
  if (j==2) {stat_lat_df<-hexp_lat_df; y_low<-0.74; y_high<-0.86; stat_name<-"Expected Heterozygosity"}
  
  colnames(stat_lat_df) <- c("Mean_Lat", stat_name)
  rownames(stat_lat_df) <- butternut_24pop_names
  stat_lat_df$Color <- NA
  stat_lat_df[1:6,3] <- "firebrick1"
  stat_lat_df[c(8,11),3] <- "lightsalmon"
  stat_lat_df[c(7,9:10),3] <- "firebrick4"
  stat_lat_df[12:24,3] <- "dodgerblue"
  
  ##matrix for all lat combos 
  r2_pvalue_df <- matrix(nrow = 4, ncol = 4)
  
  #This loop will go over the four possibilities (the list which_pops), with and without NB and Wisconsin
  #On loops 1 and 3 it will create the PDF, e.g. with and without NB 
  #On loops and 4 it prints the without Wisconsin on top of the existing plot
  for (i in 1:4){
    if (i==1&j==1) pdf("genetic_analyses_results\\all_rich_lat_quad.pdf", width = 8, height = 6)
    if (i==1&j==2) pdf("genetic_analyses_results\\hexp_lat_quad.pdf", width = 8, height = 6)
    if (i==3&j==1) pdf("genetic_analyses_results\\all_rich_lat_quad_Wo_NB.pdf", width = 8, height = 6)
    if (i==3&j==2) pdf("genetic_analyses_results\\hexp_lat_quad_Wo_NB.pdf", width = 8, height = 6)
    
    ##create plot
    if (i==1|i==3) plot(stat_lat_df[which_pops[[i]],2]~stat_lat_df[which_pops[[i]],1], 
                        col = stat_lat_df[which_pops[[i]],3], pch = 17, 
                        main = paste(stat_name,"Compared to Mean Latitude"), ylab = stat_name, 
                        xlab = "Mean Latitude", cex = (butternut_poppr[which_pops[[i]],2]/50), ylim = c(y_low,y_high))
    ##write text
    if (i==1|i==3)text(stat_lat_df[which_pops[[i]],2]~stat_lat_df[which_pops[[i]],1], 
                       labels = butternut_24pop_names[which_pops[[i]]], cex = 0.8, pos = 1)
    
    #calculate data points 
    lat_points <- stat_lat_df[which_pops[[i]],1]
    lat_points2 <- lat_points^2
    
    ###Quadratic linear model
    genstat_quad_model_lm <-lm(stat_lat_df[which_pops[[i]],2] ~ 
                                 lat_points + lat_points2)
    
    ##Calculate summary documents 
    genstat_quad_model_lm_sum <- summary(genstat_quad_model_lm)
    
    ##calculate y-axis values 
    quad_values <- seq(min(stat_lat_df[,1]), max(stat_lat_df[,1]), 0.5)
    quad_counts <- predict(genstat_quad_model_lm, list(lat_points=quad_values, lat_points2=quad_values^2))
    
    ##create r2 and p values
    genstat_r2 <- genstat_quad_model_lm_sum$adj.r.squared
    genstat_pvalue <- genstat_quad_model_lm_sum$coefficients[2,4]
    genstat_rp <- vector('expression',2)
    genstat_rp[1] <- substitute(expression(italic(R)^2 == MYVALUE), 
                                list(MYVALUE = format(genstat_r2,dig=3)))[2]
    genstat_rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                               list(MYOTHERVALUE = format(genstat_pvalue, digits = 2)))[2]
    
    r2_pvalue_df[i,3] <- genstat_r2 
    r2_pvalue_df[i,4] <- genstat_pvalue
    
    # abline(genstat_lm, col = line_colors[i])
    lines(quad_values, quad_counts, col = line_colors[i], lwd = 3)
    ##Write legend
    legend(loc_legend[i], legend = genstat_rp, bty = 'n', border = "black", 
           pt.cex = 1, cex = 0.8, pch = 17, col = line_colors[i],
           title = title_legend[i])
    legend('bottom', legend = c("New Brunswick", "Ontario", "Quebec","United States"), 
           pch = 17, col = c("firebrick1", "firebrick4","lightsalmon","dodgerblue"))
    if(i==2) dev.off()
  }
  dev.off()
}
dev.off()



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
points_counts <- predict(quadratic_model_lm, list(points=points_values, points2=points_values^2))

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

##write out data frame of statistics 
write.csv(lat_all_rich_df, "genetic_analyses_results\\r2_pvalue_lat_all_rich_df.csv")

##write out allelic richness data frame
write.csv(all_rich_lat_df, "genetic_analyses_results\\all_rich_lat_df.csv")
