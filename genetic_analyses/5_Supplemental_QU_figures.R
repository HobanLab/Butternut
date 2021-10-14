######################################
########### Load Libraries ###########
######################################

library(adegenet)
library(poppr)
library(dplyr)
library(geosphere)
library(hierfstat)

##########################################
############# Set directories ############
##########################################
##set directory to all butternut files 
butternut_drive <- "C:\\Users\\eschumacher\\Documents\\GitHub\\Butternut"

##########################################
############ Load in files ###############
##########################################
setwd(butternut_drive)

##load in 26 pop file 
butternutgen_26pop <- read.genepop("Genetic_Analyses\\data_files\\after_reorg\\butternutgen_26pop.gen", ncode = 3)

##load in lon/lat 
butternut_26pop_lonlat <- read.csv("Genetic_Analyses\\data_files\\after_reorg\\butternut_26pop_lonlat.csv")

##create population name data files 
butternut_26pop_names <- unique(butternut_26pop_lonlat$Pop)

##rename populations
levels(butternutgen_26pop@pop) <- butternut_26pop_names
##name individuals 
rownames(butternutgen_26pop@tab) <- butternut_26pop_lonlat$Ind

##load range buffer for butternut
butternut_buffer <- readOGR(dsn = paste0(butternut_drive,"\\Genetic_Analyses\\data_files\\geographic_files") , 
                            layer = "butternut_buffer")

############################################
############# Lon/Lat ######################
############################################
####Geographic analyses for each population
##calculate mean longitude and latitude for each population
butternut_mean_lon <- matrix()
butternut_mean_lat <- matrix()

##loops for mean lat/lon
for(pop in butternut_26pop_names){
  
  butternut_mean_lon[pop] <- mean(butternut_26pop_lonlat[butternut_26pop_lonlat$Pop == pop,][,3])  
  
}

for(pop in butternut_26pop_names){
  
  butternut_mean_lat[pop] <- mean(butternut_26pop_lonlat[butternut_26pop_lonlat$Pop == pop,][,4])  
  
}

##now create a matrix 
butternut_26pop_coords <- cbind(butternut_mean_lon, butternut_mean_lat)[-1,]

#####calculate distance to range edge 
#project range extent buffer to the same extent as population lon/lat
butternut_buffer_trans <- sp::spTransform(butternut_buffer, "+proj=longlat +ellps=WGS84 +datum=WGS84")

##now calculate distance to butternut range edge 
butternut_dist <- dist2Line(butternut_26pop_coords, butternut_buffer_trans)

###calculate genetic diversity stats for each population 
##calculate allelic richness and heterozygosity
butternut_26pop_allrich <- colMeans(allelic.richness(butternutgen_26pop)$Ar)

#create poppr file
butternut_26pop_poppr <- poppr(butternutgen_26pop)

##########create geographic/genetic diversity regression 
##combine stats 
butternut_26pop_gendiv_geo <- cbind(butternut_26pop_coords, data.frame(butternut_dist)[,1]/1000, butternut_26pop_allrich, 
                                    butternut_26pop_poppr[1:26,10])
##add rownames and colnames 
colnames(butternut_26pop_gendiv_geo) <- c("Mean_Lon", "Mean_Lat", "Dist_Edge", "All_Rich", "HExp")
rownames(butternut_26pop_gendiv_geo) <- butternut_26pop_names

##conver to a data frame 
butternut_26pop_gendiv_geo <- data.frame(butternut_26pop_gendiv_geo)

##add a color column 
butternut_26pop_gendiv_geo <- butternut_26pop_gendiv_geo %>% mutate(Col = NA)

butternut_26pop_gendiv_geo[1:6,6] <- "firebrick1"
butternut_26pop_gendiv_geo[7:11,6] <- "firebrick4"
butternut_26pop_gendiv_geo[12:13,6] <- "lightsalmon"
butternut_26pop_gendiv_geo[14:26,6] <- "dodgerblue"

##also create data file without WI pops 
butternut_26pop_gendiv_geo_red <- butternut_26pop_gendiv_geo[-c(18,24,25),]

########################################################
################# Run regression code ##################
########################################################
##compare genetic diversity and geographic info
linear_lat_allrich_lm <- lm(butternut_26pop_gendiv_geo$All_Rich~butternut_26pop_gendiv_geo$Mean_Lat)
linear_lat_allrich_lm_sum <- summary(linear_lat_allrich_lm)

##now compare with reduced dataset 
linear_lat_allrich_red_lm <- lm(butternut_26pop_gendiv_geo_red$All_Rich~butternut_26pop_gendiv_geo_red$Mean_Lat)
linear_lat_allrich_red_lm_sum <- summary(linear_lat_allrich_red_lm)

#allelic richness and hexp compared to latitude 
##create a linear model
#calculate data points 
lat_allrich_points <- butternut_26pop_gendiv_geo$Mean_Lat
lat_allrich_points2 <- lat_points^2

###Quadratic linear model
lat_allrich_quad_model_lm <-lm(butternut_26pop_gendiv_geo$All_Rich ~ 
                                 lat_allrich_points + lat_allrich_points2)

lat_allrich_quad_model_lm_sum <- summary(lat_allrich_quad_model_lm)

##create a data frame to store the rp values 
butternut_26pop_allrich_rp <- matrix(nrow = 4, ncol = 4)

##create regression for reduced dataset 
lat_allrich_red_points <- butternut_26pop_gendiv_geo_red$Mean_Lat
lat_allrich_red_points2 <- lat_allrich_red_points^2

###Quadratic linear model
lat_allrich_red_quad_model_lm <-lm(butternut_26pop_gendiv_geo_red$All_Rich ~ 
                                 lat_allrich_red_points + lat_allrich_red_points2)

lat_allrich_red_quad_model_lm_sum <- summary(lat_allrich_red_quad_model_lm)

##distance to edge regression 
dist_allrich_linear_lm <- lm(butternut_26pop_gendiv_geo$All_Rich~butternut_26pop_gendiv_geo$Dist_Edge)
dist_allrich_linear_lm_sum <- summary(dist_allrich_linear_lm)

##distance to edge without WI pops 
dist_allrich_linear_red_lm <- lm(butternut_26pop_gendiv_geo_red$All_Rich~butternut_26pop_gendiv_geo_red$Dist_Edge)
dist_allrich_linear_red_lm_sum <- summary(dist_allrich_linear_red_lm)

##put in table
rownames(butternut_26pop_allrich_rp) <- c("With WI pops", "Without WI pops", "With WI pops", "Without WI pops")
colnames(butternut_26pop_allrich_rp) <- c("R2","p-value","R2", "p-value")

#fill table
butternut_26pop_allrich_rp[1,1:2] <- c(as.numeric(linear_lat_allrich_lm_sum[9]),
                                       as.numeric(linear_lat_allrich_lm_sum$coefficients[2,4]))

butternut_26pop_allrich_rp[2,1:2] <- c(as.numeric(linear_lat_allrich_red_lm_sum[9]),
                                       as.numeric(linear_lat_allrich_red_lm_sum$coefficients[2,4]))

butternut_26pop_allrich_rp[1,3:4] <- c(as.numeric(lat_allrich_quad_model_lm_sum[9]), 
                                       as.numeric(lat_allrich_quad_model_lm_sum$coefficients[2,4]))

butternut_26pop_allrich_rp[2,3:4] <- c(as.numeric(lat_allrich_red_quad_model_lm_sum[9]), 
                                       as.numeric(lat_allrich_red_quad_model_lm_sum$coefficients[2,4]))

butternut_26pop_allrich_rp[3,1:2] <- c(as.numeric(dist_allrich_linear_lm_sum[9]), 
                                       as.numeric(dist_allrich_linear_lm_sum$coefficients[2,4]))

butternut_26pop_allrich_rp[4,1:2] <- c(as.numeric(dist_allrich_linear_red_lm_sum[9]), 
                                       as.numeric(dist_allrich_linear_red_lm_sum$coefficients[2,4]))

write.csv(butternut_26pop_allrich_rp, "butternut_26pop_allrich_rp.csv")

#################hexp table 
##compare genetic diversity and geographic info
linear_lat_hexp_lm <- lm(butternut_26pop_gendiv_geo$HExp~butternut_26pop_gendiv_geo$Mean_Lat)
linear_lat_hexp_lm_sum <- summary(linear_lat_hexp_lm)

##linear and lat for hexp 
linear_lat_hexp_red_lm <- lm(butternut_26pop_gendiv_geo_red$HExp~butternut_26pop_gendiv_geo_red$Mean_Lat)
linear_lat_hexp_red_lm_sum <- summary(linear_lat_hexp_red_lm)

##run quadratic regression on mean lat and hexp 
#calculate data points 
quad_lat_hexp_points <- butternut_26pop_gendiv_geo$Mean_Lat
quad_lat_hexp_points2 <- quad_lat_hexp_points^2

###Quadratic linear model
quad_lat_hexp_model_lm <-lm(butternut_26pop_gendiv_geo$HExp ~ 
                              quad_lat_hexp_points + quad_lat_hexp_points2)

##Calculate summary documents 
quad_lat_hexp_lm_sum <- summary(quad_lat_hexp_model_lm)

####calculate reduced data frame with the quadratic model 
quad_lat_hexp_red_points <- butternut_26pop_gendiv_geo_red$Mean_Lat
quad_lat_hexp_red_points2 <- quad_lat_hexp_red_points^2

##calculate reduced quadratic model 
quad_lat_hexp_model_red_lm <- lm(butternut_26pop_gendiv_geo_red$HExp ~
                                   quad_lat_hexp_red_points + quad_lat_hexp_red_points2)
##made a summary of the regression 
quad_lat_hexp_model_red_lm_sum <- summary(quad_lat_hexp_model_red_lm)

###########dist to range edge and hexp 
linear_dist_hexp_lm <- lm(butternut_26pop_gendiv_geo$HExp~butternut_26pop_gendiv_geo$Dist_Edge)
linear_dist_hexp_lm_sum <- summary(linear_dist_hexp_lm)

##now do with reduced 
linear_dist_hexp_red_lm <- lm(butternut_26pop_gendiv_geo_red$HExp~butternut_26pop_gendiv_geo_red$Dist_Edge)
linear_dist_hexp_red_lm_sum <- summary(linear_dist_hexp_red_lm)

#####create data frame 
butternut_26pop_hexp_rp_df <- matrix(nrow = 4, ncol = 4)

##make a df
butternut_26pop_hexp_rp_df[1,1:2] <- c(as.numeric(linear_lat_hexp_lm_sum[9]),
                                       as.numeric(linear_lat_hexp_lm_sum$coefficients[2,4]))

butternut_26pop_hexp_rp_df[2,1:2] <- c(as.numeric(linear_lat_hexp_red_lm_sum[9]),
                                       as.numeric(linear_lat_hexp_red_lm_sum$coefficients[2,4]))

butternut_26pop_hexp_rp_df[3,1:2] <- c(as.numeric(quad_lat_hexp_lm_sum[9]),
                                       as.numeric(quad_lat_hexp_lm_sum$coefficients[2,4]))

butternut_26pop_hexp_rp_df[4,1:2] <- c(as.numeric(quad_lat_hexp_model_red_lm_sum[9]),
                                       as.numeric(quad_lat_hexp_model_red_lm_sum$coefficients[2,4]))

butternut_26pop_hexp_rp_df[1,3:4] <- c(as.numeric(linear_dist_hexp_lm_sum[9]),
                                       as.numeric(linear_dist_hexp_lm_sum$coefficients[2,4]))

butternut_26pop_hexp_rp_df[2,3:4] <- c(as.numeric(linear_dist_hexp_red_lm_sum[9]),
                                       as.numeric(linear_dist_hexp_red_lm_sum$coefficients[2,4]))

write.csv(signif(butternut_26pop_hexp_rp_df, 3), "butternut_26pop_hexp_rp_df.csv")

##make plots 
pdf("allrich_lat.pdf", width = 8, height = 6)
plot(butternut_26pop_gendiv_geo$All_Rich~butternut_26pop_gendiv_geo$Mean_Lat, pch = 17, col = butternut_26pop_gendiv_geo$Col,
     cex = butternut_26pop_poppr[1:26,2]/20, main = "Allelic Richness Compared to Mean Population Latitude",
     xlab = "Mean Population Latitude", ylab = "Allelic Richness", cex.axis = 1.5, cex.lab = 1.5, 
     xlim = c(34,48), ylim = c(5,10))
text(butternut_26pop_gendiv_geo$All_Rich[18]~butternut_26pop_gendiv_geo$Mean_Lat[18], labels = "WI3", pos = 3, cex = 1.2)
text(butternut_26pop_gendiv_geo$All_Rich[24]~butternut_26pop_gendiv_geo$Mean_Lat[24], labels = "WI1", pos = 3, cex = 1.2)
text(butternut_26pop_gendiv_geo$All_Rich[25]~butternut_26pop_gendiv_geo$Mean_Lat[25], labels = "WI2", pos = 3, cex = 1.2)

##add legend 
legend("topleft", legend = genstat_rp, bty = 'n', border = "black", 
       pt.cex = 1.5, cex = 1.2, pch = 17, col = "darkslategray2",
       title = 'With WI populations')

legend("bottomleft", legend = genstat_red_rp, bty = 'n', border = "black", 
       pt.cex = 1.5, cex = 1.2, pch = 17, col = "darkcyan",
       title = 'Without WI populations')

legend('bottom', legend = c("New Brunswick", "Ontario", "Quebec","United States"), 
       pch = 17, col = c("firebrick1", "firebrick4","lightsalmon","dodgerblue"))

##plot quadratic line on the plot
lines(quad_values, quad_counts, lwd = 2, col = "darkslategray2")
##plot quadratic line on the plot
lines(quad_values_red, quad_counts_red, lwd = 2, col = "darkcyan")
dev.off()

###########hexp and lat 
###rp hexp 
#calculate data points 
lat_points <- butternut_26pop_gendiv_geo$Mean_Lat
lat_points2 <- lat_points^2

###Quadratic linear model
genstat_hexp_lat_quad_model_lm <-lm(butternut_26pop_gendiv_geo$HExp ~ 
                             lat_points + lat_points2)

##Calculate summary documents 
genstat_quad_hexp_lat_model_lm_sum <- summary(genstat_hexp_lat_quad_model_lm)

##calculate y-axis values 
quad_values <- seq(min(butternut_26pop_gendiv_geo[,2]), max(butternut_26pop_gendiv_geo[,2]), 0.5)
quad_counts <- predict(genstat_hexp_lat_quad_model_lm, list(lat_points=quad_values, lat_points2=quad_values^2))

##create r2 and p values
genstat_hexp_lat_r2 <- genstat_quad_hexp_lat_model_lm_sum$adj.r.squared
genstat_hexp_lat_pvalue <- genstat_quad_hexp_lat_model_lm_sum$coefficients[2,4]
genstat_hexp_lat_rp <- vector('expression',2)
genstat_rp[1] <- substitute(expression(italic(R)^2 == MYVALUE), 
                            list(MYVALUE = format(genstat_hexp_lat_r2,dig=3)))[2]
genstat_rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                           list(MYOTHERVALUE = format(genstat_hexp_lat_pvalue, digits = 2)))[2]

##without WI 
#calculate data points 
hexp_lat_points_red <- butternut_26pop_gendiv_geo_red$Mean_Lat
hexp_lat_points2_red <- hexp_lat_points_red^2

###Quadratic linear model
genstat_hexp_lat_quad_model_red_lm <-lm(butternut_26pop_gendiv_geo_red$HExp ~ 
                                 hexp_lat_points_red + hexp_lat_points2_red)

##Calculate summary documents 
genstat_hexp_lat_quad_model_red_lm_sum <- summary(genstat_hexp_lat_quad_model_red_lm)

##calculate y-axis values 
quad_hexp_lat_values_red <- seq(min(butternut_26pop_gendiv_geo_red[,2]), max(butternut_26pop_gendiv_geo_red[,2]), 0.5)
quad_hexp_lat_counts_red <- predict(genstat_hexp_lat_quad_model_red_lm, list(hexp_lat_points_red=quad_hexp_lat_values_red, 
                                                                             hexp_lat_points2_red=quad_hexp_lat_values_red^2))

##create r2 and p values
genstat_hexp_lat_red_r2 <- genstat_hexp_lat_quad_model_red_lm_sum$adj.r.squared
genstat_hexp_lat_red_pvalue <- genstat_hexp_lat_quad_model_red_lm_sum$coefficients[2,4]
genstat_hexp_lat_red_rp <- vector('expression',2)
genstat_hexp_lat_red_rp[1] <- substitute(expression(italic(R)^2 == MYVALUE), 
                                list(MYVALUE = format(genstat_hexp_lat_red_r2,dig=3)))[2]
genstat_hexp_lat_red_rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                               list(MYOTHERVALUE = format(genstat_hexp_lat_red_pvalue, digits = 2)))[2]

##make plots 
pdf("hexp_lat.pdf", width = 8, height = 6)
plot(butternut_26pop_gendiv_geo$HExp~butternut_26pop_gendiv_geo$Mean_Lat, pch = 17, col = butternut_26pop_gendiv_geo$Col,
     cex = butternut_26pop_poppr[1:26,2]/20, main = "Expected Heterozygosity Compared to Mean Population Latitude",
     xlab = "Mean Population Latitude", ylab = "Expected Heterozygosity", cex.axis = 1.5, cex.lab = 1.5, 
     xlim = c(34,48), ylim = c(0.71,0.9))
text(butternut_26pop_gendiv_geo$HExp[18]~butternut_26pop_gendiv_geo$Mean_Lat[18], labels = "WI3", pos = 3, cex = 1.2)
text(butternut_26pop_gendiv_geo$HExp[24]~butternut_26pop_gendiv_geo$Mean_Lat[24], labels = "WI1", pos = 3, cex = 1.2)
text(butternut_26pop_gendiv_geo$HExp[25]~butternut_26pop_gendiv_geo$Mean_Lat[25], labels = "WI2", pos = 3, cex = 1.2)

##add legend 
legend("topleft", legend = genstat_rp, bty = 'n', border = "black", 
       pt.cex = 1.5, cex = 1.2, pch = 17, col = "darkslategray2",
       title = 'With WI populations')

legend("bottomleft", legend = genstat_hexp_lat_red_rp, bty = 'n', border = "black", 
       pt.cex = 1.5, cex = 1.2, pch = 17, col = "darkcyan",
       title = 'Without WI populations')

legend('bottom', legend = c("New Brunswick", "Ontario", "Quebec","United States"), 
       pch = 17, col = c("firebrick1", "firebrick4","lightsalmon","dodgerblue"))

##plot quadratic line on the plot
lines(quad_values, quad_counts, lwd = 2, col = "darkslategray2")
##plot quadratic line on the plot
lines(quad_hexp_lat_values_red, quad_hexp_lat_counts_red, lwd = 2, col = "darkcyan")
dev.off()

######################distance to range edge figures
###linear regression 
dist_allrich_lm <- lm(butternut_26pop_gendiv_geo$All_Rich ~ butternut_26pop_gendiv_geo$Dist_Edge)
dist_allrich_lm_sum <- summary(dist_allrich_lm)

##linear regression without WI 
dist_allrich_red_lm <- lm(butternut_26pop_gendiv_geo_red$All_Rich ~ butternut_26pop_gendiv_geo_red$Dist_Edge)
dist_allrich_red_lm_sum <- summary(dist_allrich_red_lm)

##create r2 and p values
dist_allrich_r2 <- dist_allrich_lm_sum$adj.r.squared
dist_allrich_pvalue <- dist_allrich_lm_sum$coefficients[2,4]
dist_allrich_rp <- vector('expression',2)
dist_allrich_rp[1] <- substitute(expression(italic(R)^2 == MYVALUE), 
                                list(MYVALUE = format(dist_allrich_r2,dig=3)))[2]
dist_allrich_rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                               list(MYOTHERVALUE = format(dist_allrich_pvalue, digits = 2)))[2]

##create r2 and pvalues 
dist_allrich_red_r2 <- dist_allrich_red_lm_sum$adj.r.squared
dist_allrich_red_pvalue <- dist_allrich_red_lm_sum$coefficients[2,4]
dist_allrich_red_rp <- vector('expression',2)
dist_allrich_red_rp[1] <- substitute(expression(italic(R)^2 == MYVALUE), 
                                 list(MYVALUE = format(dist_allrich_red_r2,dig=3)))[2]
dist_allrich_red_rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                                list(MYOTHERVALUE = format(dist_allrich_pvalue, digits = 2)))[2]

##make plots 
pdf("allrich_distedge.pdf", width = 8, height = 6)
plot(butternut_26pop_gendiv_geo$All_Rich~butternut_26pop_gendiv_geo$Dist_Edge, pch = 17, col = butternut_26pop_gendiv_geo$Col,
     cex = butternut_26pop_poppr[1:26,2]/20, main = "Allelic Richness Compared to Distance to Range Edge (km)",
     xlab = "Distance to Range Edge (km)", ylab = "Allelic Richness", cex.axis = 1.5, cex.lab = 1.5, 
     xlim = c(0,600), ylim = c(5,10))
text(butternut_26pop_gendiv_geo$All_Rich[18]~butternut_26pop_gendiv_geo$Dist_Edge[18], labels = "WI3", pos = 3, cex = 1.2)
text(butternut_26pop_gendiv_geo$All_Rich[24]~butternut_26pop_gendiv_geo$Dist_Edge[24], labels = "WI1", pos = 3, cex = 1.2)
text(butternut_26pop_gendiv_geo$All_Rich[25]~butternut_26pop_gendiv_geo$Dist_Edge[25], labels = "WI2", pos = 3, cex = 1.2)

abline(dist_allrich_lm, col = "dodgerblue4", lwd = 2)
abline(dist_allrich_red_lm, col = "darkorchid", lwd = 2)

##add legend 
legend("topleft", legend = dist_allrich_rp, bty = 'n', border = "black", 
       pt.cex = 1.5, cex = 1.2, pch = 17, col = "dodgerblue4",
       title = 'With WI populations')

legend("bottomleft", legend = genstat_hexp_lat_red_rp, bty = 'n', border = "black", 
       pt.cex = 1.5, cex = 1.2, pch = 17, col = "darkorchid",
       title = 'Without WI populations')

legend('bottom', legend = c("New Brunswick", "Ontario", "Quebec","United States"), 
       pch = 17, col = c("firebrick1", "firebrick4","lightsalmon","dodgerblue"))

dev.off()

#################dist to range edge hexp
##create linear regressions 
dist_hexp_lm <- lm(butternut_26pop_gendiv_geo$HExp ~ butternut_26pop_gendiv_geo$Dist_Edge)
dist_hexp_lm_sum <- summary(dist_hexp_lm)

##linear regression without WI 
dist_hexp_red_lm <- lm(butternut_26pop_gendiv_geo_red$HExp ~ butternut_26pop_gendiv_geo_red$Dist_Edge)
dist_hexp_red_lm_sum <- summary(dist_hexp_red_lm)

##create r2 and p values
dist_hexp_r2 <- dist_hexp_lm_sum$adj.r.squared
dist_hexp_pvalue <- dist_hexp_lm_sum$coefficients[2,4]
dist_hexp_rp <- vector('expression',2)
dist_hexp_rp[1] <- substitute(expression(italic(R)^2 == MYVALUE), 
                                 list(MYVALUE = format(dist_hexp_r2,dig=3)))[2]
dist_hexp_rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                                list(MYOTHERVALUE = format(dist_hexp_pvalue, digits = 2)))[2]

##create r2 and pvalues 
dist_hexp_red_r2 <- dist_hexp_red_lm_sum$adj.r.squared
dist_hexp_red_pvalue <- dist_hexp_red_lm_sum$coefficients[2,4]
dist_hexp_red_rp <- vector('expression',2)
dist_hexp_red_rp[1] <- substitute(expression(italic(R)^2 == MYVALUE), 
                                     list(MYVALUE = format(dist_hexp_red_r2,dig=3)))[2]
dist_hexp_red_rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                                    list(MYOTHERVALUE = format(dist_hexp_red_pvalue, digits = 2)))[2]

##make plots 
pdf("hexp_distedge.pdf", width = 8, height = 6)
plot(butternut_26pop_gendiv_geo$HExp~butternut_26pop_gendiv_geo$Dist_Edge, pch = 17, col = butternut_26pop_gendiv_geo$Col,
     cex = butternut_26pop_poppr[1:26,2]/20, main = "Expected Heterozygosity Compared to Distance to Range Edge (km)",
     xlab = "Distance to Range Edge (km)", ylab = "Expected Heterozygosity", cex.axis = 1.5, cex.lab = 1.5, 
     xlim = c(0,600), ylim = c(0.71, 0.9))
text(butternut_26pop_gendiv_geo$HExp[18]~butternut_26pop_gendiv_geo$Dist_Edge[18], labels = "WI3", pos = 3, cex = 1.2)
text(butternut_26pop_gendiv_geo$HExp[24]~butternut_26pop_gendiv_geo$Dist_Edge[24], labels = "WI1", pos = 3, cex = 1.2)
text(butternut_26pop_gendiv_geo$HExp[25]~butternut_26pop_gendiv_geo$Dist_Edge[25], labels = "WI2", pos = 3, cex = 1.2)

abline(dist_hexp_lm, col = "dodgerblue4", lwd = 2)
abline(dist_hexp_red_lm, col = "darkorchid", lwd = 2)

##add legend 
legend("topleft", legend = dist_hexp_rp, bty = 'n', border = "black", 
       pt.cex = 1.5, cex = 1.2, pch = 17, col = "dodgerblue4",
       title = 'With WI populations')

legend("bottomleft", legend = dist_hexp_red_rp, bty = 'n', border = "black", 
       pt.cex = 1.5, cex = 1.2, pch = 17, col = "darkorchid",
       title = 'Without WI populations')

legend('bottom', legend = c("New Brunswick", "Ontario", "Quebec","United States"), 
       pch = 17, col = c("firebrick1", "firebrick4","lightsalmon","dodgerblue"))

dev.off()

########################################
############## PCOA ####################
########################################
##first, code to reduce by individuals 
butternut_26pop_relate_df <- Demerelate(butternut_26pop_lonlat[,-c(3:4)], object = T, value = "loiselle")

##now identify how many individuals have greater than 25% relatedness = half siblings
butternut_halfsib_names <- names(which(unlist(butternut_26pop_relate_df$Empirical_Relatedness) > 0.25))

##then use this to create a document which has all the unique individual numbers for every highly related individuals
butternut_halfsib_names_cleanfront <- gsub("^.*\\.","", butternut_halfsib_names)

butternut_halfsib_names_cleanback <- gsub("^.*\\_","", butternut_halfsib_names_cleanfront)

relate_ind_remove <- unique(butternut_halfsib_names_cleanback)

##then subset genind file
butternutgen_26pop_relatedness_reduced <- butternutgen_26pop[!rownames(butternutgen_26pop@tab) %in% relate_ind_remove,]

##subset data frame 
butternut_26pop_relatedness_reduced_df <- butternut_26pop_lonlat[!butternut_26pop_lonlat$Ind %in% relate_ind_remove,]

##name pops in the genind file 
levels(butternutgen_26pop_relatedness_reduced@pop) <- butternut_26pop_names

##export genind 2 genalex 
butternut_26pop_red_genalex <- genind2genalex(butternutgen_26pop_relatedness_reduced, 
                                              filename = 'butternutgen_26pop_relatedness_reduced.csv')

###########make PCoA 
##run PCOA analysis 
##now reorg into a genepop file 
butternut_26pop_relate <- genind2genpop(butternutgen_26pop_relatedness_reduced)

##run the PCA function
butternut_reorg_pco <- dudi.pco(dist.genpop(butternut_26pop_relate, meth = 2), nf = 2, scannf = FALSE)

##format df
rownames(butternut_reorg_pco$li) <- butternut_26pop_names

##re-org data for color coding
butternut_pco_nb <- rbind(butternut_reorg_pco$li[c(1:6),1:2])
butternut_pco_ot <- rbind(butternut_reorg_pco$li[c(7:11),1:2])
butternut_pco_qu <- rbind(butternut_reorg_pco$li[c(12:13),1:2])
butternut_pco_us <- rbind(butternut_reorg_pco$li[c(14:26),1:2])

##Calculate percent variation explained 
sum_eig <- sum(butternut_reorg_pco$eig)
pc1 <- (butternut_reorg_pco$eig[[1]]/sum_eig)*100
pc2 <- (butternut_reorg_pco$eig[[2]]/sum_eig)*100

##PCOA of the reorg
pdf("Genetic_Analyses\\genetic_analyses_results\\PCoA.pdf", width = 8, height = 6)
##plot New Brunswick populations
plot(butternut_pco_nb$A1, butternut_pco_nb$A2, pch = 17, 
     xlab = paste0("PC1", sep = " ", "(",round(pc1, digits = 1), "%", ")"), 
     ylab = paste0("PC2", sep = " ", "(",round(pc2, digits = 1), "%",")"), 
     col = "firebrick1", xlim = c(-0.4, 0.4), 
     ylim = c(-0.4, 0.4), cex = 1.2)
##make axes
abline(h = 0)
abline(v = 0)
##plot Ontario populations
points(butternut_pco_ot$A1, butternut_pco_ot$A2, pch = 17, col = "firebrick4", cex = 1.2)

##plot quebec populations 
points(butternut_pco_qu$A1, butternut_pco_qu$A2, pch = 17, col = "lightsalmon", cex = 1.2)

##plot US populations 
points(butternut_pco_us$A1, butternut_pco_us$A2, col = "dodgerblue", pch = 17, cex = 1.2)

##label Wisconsin populations 
text(butternut_reorg_pco$li$A1[24], butternut_reorg_pco$li$A2[24], labels = "WI1", pos = 3, cex = 1.2)
text(butternut_reorg_pco$li$A1[25], butternut_reorg_pco$li$A2[25], labels = "WI2", pos = 1, cex = 1.2)
text(butternut_reorg_pco$li$A1[18], butternut_reorg_pco$li$A2[18], labels = "WI3", pos = 3, cex = 1.2)


##create legend 
legend('topleft', pch = 17, col = c("firebrick1","firebrick4", "lightsalmon", "dodgerblue"), 
       legend = c("New Brunswick", "Ontario", "Quebec", "United States"), cex = 1.2)

dev.off()

##run session info 
sessionInfo()



