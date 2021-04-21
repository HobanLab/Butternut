##########################
######## Libraries #######
##########################

library(geosphere)
library(raster)

#####################################
############ Directories ############
#####################################
butternut_drive <- "C:\\Users\\eschumacher\\Documents\\GitHub\\butternut"

##############################################################
####################### Load Files  ##########################
##############################################################
setwd(butternut_drive)

##load in butternut buffer shapefile 
butternut_buffer <- shapefile("data_files\\geographic_files\\butternut_buffer")

##load mean longitude and latitude document for all populations 
butternut_coord_df <- read.csv("data_files\\geographic_files\\butternut_coord_df.csv")

##load relatedness file to name individuals
reorg_relatedness <- read.csv("data_files\\after_reorg\\reorg_relatedness.csv")

####################################################
############# Calculate distance to edge ###########
####################################################
##prep data frame for spatial analysis
##create list of the 24 populations 
butternut_24pop_names <- unique(reorg_relatedness$Pop)

##get butternut 24 population names 
butternut_coord_df[,1] <- butternut_24pop_names

##name columns 
colnames(butternut_coord_df) <- c("Pop","Longitude","Latitude","Col")

##make spatial
coordinates(butternut_coord_df) <- c('Longitude', 'Latitude')
proj4string(butternut_coord_df) <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84")

##change proj4string 
butternut_buffer_WGS84 <- spTransform(butternut_buffer, CRS("+proj=longlat +ellps=WGS84 +datum=WGS84"))

#####Calculate distance to edge in meters 
dist_to_edge <- dist2Line(butternut_coord_df, butternut_buffer_WGS84, distfun = distGeo)

##put butternut coord df back into a data frame 
butternut_lonlat_df <- data.frame(butternut_coord_df)

####Now combine into a data frame, divide by 1000 to get into km 
butternut_dist_df <- data.frame(cbind(butternut_lonlat_df[,1:3], (dist_to_edge[,1]/1000)))

##add column names 
colnames(butternut_dist_df) <- c("Pop","Mean_Lon", "Mean_Lat", "Dist_To_Edge")

##add a column for colors by location 
butternut_dist_df$Col <- "NA"

##add population colors for graphing 
butternut_dist_df[1:6,5] <- "firebrick1"
butternut_dist_df[c(8,11),5] <- "lightsalmon"
butternut_dist_df[c(7,9:10),5] <- "firebrick4"
butternut_dist_df[12:24,5] <- "dodgerblue"

##write out distance to edge csv
write.csv(butternut_dist_df, "data_files\\geographic_files\\butternut_dist_edge_df.csv")
