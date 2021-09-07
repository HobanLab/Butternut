#This code runs all of the regression analyses on geographic 
#location with butternut population's genetic diversity 

##########################
######## Libraries #######
##########################

library(adegenet)
library(poppr)
library(hierfstat)

#####################################
############ Directories ############
#####################################
butternut_drive <- "C:\\Users\\eschumacher\\Documents\\GitHub\\Butternut"

#####################################
############# Load Files ############
#####################################
setwd(butternut_drive)

##load current working reorganized genepop file 
butternutgen_reorg <- read.genepop("Genetic_Analyses\\data_files\\after_reorg\\reorg_gen_24pop.gen", ncode = 3)

##load relatedness document to name individuals 
reorg_relatedness <- read.csv("Genetic_Analyses\\data_files\\after_reorg\\reorg_relatedness.csv")

##name individuals in the genind file  
rownames(butternutgen_reorg@tab) <- reorg_relatedness$Ind

##load in mean longitude and latitude by population 
butternut_mean_lonlat <- read.csv("Genetic_Analyses\\data_files\\geographic_files\\butternut_coord_df.csv")

##get pop names 
butternut_24pop_names <- unique(reorg_relatedness$Pop)

##name butternut lon and lat 
butternut_mean_lonlat[,1] <- butternut_24pop_names

##name populations in genind file 
levels(butternutgen_reorg@pop) <- butternut_24pop_names

##create butternut poppr
butternut_poppr <- poppr(butternutgen_reorg)

###################################################################################
#### Linear Modeling on Allelic_Richness and Expected Heterozygosity ##############
###################################################################################

#We have a nested loop below.  The outer loop goes over the statistics, alleles and heterozygosity
#The inner loop goes over the four subsets of population groupings
#These are the subsets: all pops, all minus WI, all minus NB, all minus (NB+WI)
#Also includes regression type - quadratic vs. linear
#And models distance to range edge to these two statistics
which_pops<-list(1:24,c(1:15,17:21,24),7:24,c(7:15,17:21,24),
                 1:24,c(1:15,17:21,24),7:24,c(7:15,17:21,24),
                 1:24,c(1:15,17:21,24),7:24,c(7:15,17:21,24))
#The following are options used in the loop of 4
line_colors<-c("dodgerblue4","darkorchid4","dodgerblue4","darkorchid4",
               "darkslategray3", "darkseagreen4", "darkslategray3", "darkseagreen4")
loc_legend<-c("topleft","bottomleft","topleft","bottomleft", 
              "topleft","bottomleft","topleft","bottomleft",
              "topleft","bottomleft","topleft","bottomleft")
title_legend<-c("With WI Populations", "Without WI Populations","With WI Populations", "Without WI Populations",
                "With WI Populations", "Without WI Populations","With WI Populations", "Without WI Populations")

##calculate allelic richness
reorg_allrich <- colMeans(allelic.richness(butternutgen_reorg)$Ar)

##create df with allelic richness and latitude
all_rich_lat_df <- data.frame(butternut_mean_lonlat[,3], reorg_allrich)
##create data frame with heterozygosity and latitude
hexp_lat_df <- read.csv("Genetic_Analyses\\genetic_analyses_results\\hexp_lat_df.csv")[,c(2,4)]

##distance to edge data frames with genetic diversity stats
all_rich_dist_edge_df <- read.csv("Genetic_Analyses\\genetic_analyses_results\\allrich_dist_edge_df.csv")[,c(2:4)]
hexp_dist_edge_df <- read.csv("Genetic_Analyses\\genetic_analyses_results\\hexp_dist_edge_df.csv")[,c(2:4)]

##data frames to save all rp values 
allrich_rp <- matrix(nrow = 8, ncol = 4)

##create data frame to save hexp 
hexp_rp <- matrix(nrow = 8, ncol = 4)

#The loop will go over the two statistics, 1 = alleles, 2 = heterozygosity; 1 and 2 are compared to latitude
##3 and 4 are compared to distance to edge
for (j in 1:4){

    #This will set y limits and axis titles to match the statistic
    if (j==1) {stat_lat_df<-all_rich_lat_df; y_low<-5; y_high<-10; stat_name<-"Allelic Richness"}
    if (j==2) {stat_lat_df<-hexp_lat_df; y_low<-0.74; y_high<-0.86; stat_name<-"Expected Heterozygosity"}
    if (j==3) {stat_dist_df<-all_rich_dist_edge_df; y_low<-5; y_high<-10; stat_name<-"Allelic Richness"}
    if (j==4) {stat_dist_df<-hexp_dist_edge_df; y_low<-0.74; y_high<-0.86; stat_name<-"Expected Heterozygosity"}
  
    ##creates color column and names rows and columns
    colnames(stat_lat_df) <- c("Mean_Lat", stat_name)
    rownames(stat_lat_df) <- butternut_24pop_names
    stat_lat_df$Color <- NA
    stat_lat_df[1:6,3] <- "firebrick1"
    stat_lat_df[c(8,11),3] <- "lightsalmon"
    stat_lat_df[c(7,9:10),3] <- "firebrick4"
    stat_lat_df[12:24,3] <- "dodgerblue"

    #This loop will go over the eight possibilities (the list which_pops), with and without NB and Wisconsin
    #1 and 3 are linear models with latitude
    #5 and 7 are quadratic models with latitude and genetic diversity stats
    #On loops 2, 4, 6, 8 it prints without Wisconsin regression on top of the existing plot
    
    for (i in 1:8){
	    #linear regressions 
      if (i==1&j==1) pdf("Genetic_Analyses\\genetic_analyses_results\\all_rich_lat_linear.pdf", width = 8, height = 6)
	    if (i==1&j==2) pdf("Genetic_Analyses\\genetic_analyses_results\\hexp_lat_linear.pdf", width = 8, height = 6)
	    if (i==3&j==1) pdf("Genetic_Analyses\\genetic_analyses_results\\all_rich_lat_linear_Wo_NB.pdf", width = 8, height = 6)
	    if (i==3&j==2) pdf("Genetic_Analyses\\genetic_analyses_results\\hexp_lat_linear_Wo_NB.pdf", width = 8, height = 6)
      ##quadratic regressions 
      if (i==5&j==1) pdf("Genetic_Analyses\\genetic_analyses_results\\all_rich_lat_quad.pdf", width = 8, height = 6)
      if (i==5&j==2) pdf("Genetic_Analyses\\genetic_analyses_results\\hexp_lat_quad.pdf", width = 8, height = 6)
      if (i==7&j==1) pdf("Genetic_Analyses\\genetic_analyses_results\\all_rich_lat_quad_Wo_NB.pdf", width = 8, height = 6)
      if (i==7&j==2) pdf("Genetic_Analyses\\genetic_analyses_results\\hexp_lat_quad_Wo_NB.pdf", width = 8, height = 6)
      ##distance to range edge linear regressions
      if (i==1&j==3) pdf("Genetic_Analyses\\genetic_analyses_results\\all_rich_dist_edge_linear.pdf", width = 8, height = 6)
      if (i==1&j==4) pdf("Genetic_Analyses\\genetic_analyses_results\\hexp_edge_linear.pdf", width = 8, height = 6)
      if (i==3&j==3) pdf("Genetic_Analyses\\genetic_analyses_results\\all_rich_dist_edge_linear_Wo_NB.pdf", width = 8, height = 6)
      if (i==3&j==4) pdf("Genetic_Analyses\\genetic_analyses_results\\hexp_edge_linear_Wo_NB.pdf", width = 8, height = 6)

  if(j==1|j==2){
	##create plot
	if (i==1|i==3|i==5|i==7) plot(stat_lat_df[which_pops[[i]],2]~stat_lat_df[which_pops[[i]],1], 
	                              col = stat_lat_df[which_pops[[i]],3], 
	                              pch = 17, 
		 main = paste(stat_name,"Compared to Mean Latitude"), ylab = stat_name, 
		 xlab = "Mean Latitude", cex = (butternut_poppr[which_pops[[i]],2]/50), ylim = c(y_low,y_high))
  
	##write text
	if (i==1|i==3|i==5|i==7) text(stat_lat_df[which_pops[[i]],2]~stat_lat_df[which_pops[[i]],1], 
	                          labels = butternut_24pop_names[which_pops[[i]]], cex = 0.8, pos = 1)
  if(i==1|i==2|i==3|i==4){
	##create a linear model
	genstat_lm <- lm(stat_lat_df[which_pops[[i]],2]~stat_lat_df[which_pops[[i]],1])
	genstat_lm_sum <- summary(genstat_lm)

	##create r2 and p values
	genstat_r2 <- genstat_lm_sum$adj.r.squared
	genstat_pvalue <- genstat_lm_sum$coefficients[2,4]
	genstat_rp <- vector('expression',2)
	genstat_rp[1] <- substitute(expression(italic(R)^2 == MYVALUE), 
									list(MYVALUE = format(genstat_r2,dig=3)))[2]
	genstat_rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
								   list(MYOTHERVALUE = format(genstat_pvalue, digits = 2)))[2]
	
	##plot linear regression 
	abline(genstat_lm, col = line_colors[i], lwd = 2)
	
	###save rp values in a data frame 
	if(j==1){allrich_rp[i,c(1:2)] <- c(genstat_r2, genstat_pvalue)}
	if(j==2){hexp_rp[i,c(1:2)] <- c(genstat_r2, genstat_pvalue)}
	
	##generate quadratic regressions 
  }else{
  
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
    ##plot quadratic line on the plot
    lines(quad_values, quad_counts, col = line_colors[i], lwd = 2)
    
    ##save values in data frames
    if(j==1){allrich_rp[i-4,c(3:4)] <- c(genstat_r2, genstat_pvalue)}
    if(j==2){hexp_rp[i-4,c(3:4)] <- c(genstat_r2, genstat_pvalue)}
    
  }
	
	  ##Write legend
	  legend(loc_legend[i], legend = genstat_rp, bty = 'n', border = "black", 
		   pt.cex = 1, cex = 0.8, pch = 17, col = line_colors[i],
		   title = title_legend[i])
	  legend('bottom', legend = c("New Brunswick", "Ontario", "Quebec","United States"), 
		   pch = 17, col = c("firebrick1", "firebrick4","lightsalmon","dodgerblue"))
	#if(i==2|i==6) dev.off()
	  
  }else{
    
    if (i==1|i==3)
	  plot(stat_dist_df[which_pops[[i]],2]~stat_dist_df[which_pops[[i]],1], col = stat_dist_df[which_pops[[i]],3], 
	       pch = 17, 
	       main = paste(stat_name,"Compared to Distance to Edge (km)"), ylab = stat_name, 
	       xlab = "Distance to Edge (km)", cex = (butternut_poppr[which_pops[[i]],2]/50), ylim = c(y_low,y_high))
	  
	  ##write text
	  if (i==1|i==3) text(stat_dist_df[which_pops[[i]],2]~stat_dist_df[which_pops[[i]],1], 
	                      labels = butternut_24pop_names[which_pops[[i]]], cex = 0.8, pos = 1)
    
    if(i==1|i==2|i==3|i==4){
      ##create a linear model
      gendist_lm <- lm(stat_dist_df[which_pops[[i]],2]~stat_dist_df[which_pops[[i]],1])
      gendist_lm_sum <- summary(gendist_lm)
      
      ##create r2 and p values
      gendist_r2 <- gendist_lm_sum$adj.r.squared
      gendist_pvalue <- gendist_lm_sum$coefficients[2,4]
      gendist_rp <- vector('expression',2)
      gendist_rp[1] <- substitute(expression(italic(R)^2 == MYVALUE), 
                                  list(MYVALUE = format(gendist_r2,dig=3)))[2]
      gendist_rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                                 list(MYOTHERVALUE = format(gendist_pvalue, digits = 2)))[2]
      
      ##save values in data frames
      if(j==3){allrich_rp[i+4,c(1:2)] <- c(gendist_r2, gendist_pvalue)}
      if(j==4){hexp_rp[i+4,c(1:2)] <- c(gendist_r2, gendist_pvalue)}
      
      abline(gendist_lm, col = line_colors[i], lwd = 2)
	    
	    ##Write legend
	    legend(loc_legend[i], legend = gendist_rp, bty = 'n', border = "black", 
	           pt.cex = 1, cex = 0.8, pch = 17, col = line_colors[i],
	           title = title_legend[i])
	    legend('bottom', legend = c("New Brunswick", "Ontario", "Quebec","United States"), 
	           pch = 17, col = c("firebrick1", "firebrick4","lightsalmon","dodgerblue"))
	    
	        }
        }
    
      }
      
    }
   





dev.off()
dev.off()
dev.off()
dev.off()
dev.off()
dev.off()
dev.off()
dev.off()
dev.off()
dev.off()
dev.off()
dev.off()
dev.off()


##write out tables of r2 and p values 
##name rows and columns 
rownames(allrich_rp) <- c("Lat_AllPops","Lat_woWI", "Lat_WithWI_woNB","Lat_woWI_woNB",
                              "Dist_Edge_AllPops","Dist_Edge_woWI", "Dist_Edge_WithWI_woNB","Dist_Edge_woWI_woNB" )
colnames(allrich_rp) <- c("R2 Linear", "p-value Linear", "R2 Quadratic", "p-value Quadratic")

##rename columns and rownames for hexp r2 and pvalues
rownames(hexp_rp) <- c("Lat_AllPops","Lat_woWI", "Lat_WithWI_woNB","Lat_woWI_woNB",
                              "Dist_Edge_AllPops","Dist_Edge_woWI", "Dist_Edge_WithWI_woNB","Dist_Edge_woWI_woNB" )
colnames(hexp_rp) <- c("R2 Linear", "p-value Linear", "R2 Quadratic", "p-value Quadratic")

##round of significant tables
allrich_rp <- signif(allrich_rp,2)
hexp_rp <- signif(hexp_rp,2)
##write out csv
write.csv(allrich_rp, "Genetic_Analyses\\genetic_analyses_results\\allrich_rp_df.csv")
write.csv(hexp_rp, "Genetic_Analyses\\genetic_analyses_results\\hexp_rp_df.csv")

##run sessionInfo
sessionInfo()
