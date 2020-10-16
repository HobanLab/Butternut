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
dist_edge_path <- "G:\\Shared drives\\Emily_Schumacher\\butternut_publication_figures\\dist_edge\\"
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
butternut_dist <- cbind(butternut_mean_lon, butternut_mean_lat, (dist_to_edge[,1]/1000))
colnames(butternut_dist) <- c("MeanLong", "MeanLat", "Dist_To_Edge")
rownames(butternut_dist) <- butternut_24pop_names



