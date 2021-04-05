##########################
######## Libraries #######
##########################

library(diveRsity)
library(adegenet)
library(stringr)
library(tidyr)
library(hierfstat)
library(poppr)
library(Demerelate)
library(rworldmap)
library(data.table)
library(ggplot2)
library(ggrepel)
library(geosphere)
library(plotrix)
library(ggpmisc)
library(factoextra)
library(GISTools)
library(raster)
library(rgdal)
library(sp)

#####################################
############ Directories ############
#####################################

butternut_drive <- "G:\\My Drive\\Hoban_Lab_Docs\\Projects\\Butternut_JUCI\\"

##load butternut
setwd(dist_edge_path)
butternut_buffer <- shapefile("butternut_buffer")

##############################################################
#################### Load Genetic Files  #####################
##############################################################
setwd(butternut_drive)

##load current working reorganized genepop file 
butternutgen_reorg <- read.genepop("DataFiles\\24Populations\\reorg\\reorg_gen_24pop.gen", ncode = 3)

##Lon lat files
butternut_reorg_lonlat <- read.csv("DataFiles\\24Populations\\reorg\\reorg_lon_lat.csv")

###Name the reorg file
rownames(butternutgen_reorg@tab) <- butternut_reorg_lonlat$Ind

##GESTE Results
reorg_geste_fst <- read.csv(paste0(butternut_drive,"\\DataFiles\\24Populations\\reorg\\GESTE_fst.csv"))

##create population name file 
butternut_24pop_names <- unique(butternut_reorg_lonlat$Pop)

##Name levels
levels(butternutgen_reorg@pop) <- butternut_24pop_names

##calculate poppr
butternut_24pop_poppr <- poppr(butternutgen_reorg)

####################################################
########### Genetic Diversity Calcs ################
####################################################
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

###make into one data frame
butternut_coord_df <- data.frame(butternut_mean_lon, butternut_mean_lat)
rownames(butternut_coord_df) <- butternut_24pop_names
colnames(butternut_coord_df) <- c("Longitude", "Latitude")

##make spatial
coordinates(butternut_coord_df) <- c('Longitude', 'Latitude')
proj4string(butternut_coord_df) <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84")

##change proj4string 
butternut_buffer_WGS84 <- spTransform(butternut_buffer, CRS("+proj=longlat +ellps=WGS84 +datum=WGS84"))

#####Calculate distance to edge
dist_to_edge <- dist2Line(butternut_coord_df, butternut_buffer_WGS84, distfun = distGeo)

####Now combine into a data frame
butternut_dist <- data.frame(cbind(butternut_mean_lon, butternut_mean_lat, (dist_to_edge[,1]/1000)))
colnames(butternut_dist) <- c("MeanLong", "MeanLat", "Dist_To_Edge")
rownames(butternut_dist) <- butternut_24pop_names
butternut_dist$col <- NA
butternut_dist[1:6,4] <- "firebrick1"
butternut_dist[7:11,4] <- "firebrick4"
butternut_dist[12:24,4] <- "dodgerblue"

##write out distance to edge csv
write.csv(butternut_dist, paste0(dist_edge_path, "dist_edge_df.csv"))

##calculate regression
edgedist_fst_lm <- lm(reorg_geste_fst[,2]~butternut_dist[,3])
edgedist_fst_lm_sum <- summary(edgedist_fst_lm)

##create r2 and p values
edgefst_r2 <- edgedist_fst_lm_sum$adj.r.squared
edgefst_coef <- edgedist_fst_lm_sum$coefficients
edgefst_pvalue <- edgefst_coef[2,4]
edgefst_rp <- vector('expression',2)
edgefst_rp[1] <- substitute(expression(italic(R)^2 == MYVALUE), 
                           list(MYVALUE = format(edgefst_r2,dig=3)))[2]
edgefst_rp[2] <- substitute(expression(italic(p) == MYOTHERVALUE), 
                           list(MYOTHERVALUE = format(edgefst_pvalue, digits = 2)))[2]


##plot comparison 
pdf(paste0(dist_edge_path,"fst_dist_edge.pdf"), width = 8, height = 6)

plotCI(butternut_dist[,3], reorg_geste_fst[,2], ui = reorg_geste_fst[,5], li = reorg_geste_fst[,4], 
       ylim = c(0,0.15), pch = 17, xlab = "Mean Distance to Edge (km)",
       xlim = c(0,600),
       ylab = "GESTE Fst", col = butternut_dist[,4], cex = (butternut_24pop_poppr[1:24,2]/50), 
       main = "GESTE Fst Comapred to Distance to Range Edge (km)")

text(butternut_dist[,3], reorg_geste_fst[,2], labels = butternut_24pop_names, 
     cex = 0.8, pos = 3)

abline(edgedist_fst_lm, col = "dodgerblue4")

legend('topright', legend = c("New Brunswick", "Ontario", "United States"), pch = 17, col = c("firebrick1", "firebrick4","dodgerblue"))

legend('topleft', legend = edgefst_rp, bty = 'n', border = "black", pt.cex = 1, cex = 0.8)

dev.off()


