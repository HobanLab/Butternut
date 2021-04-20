##########################
######## Libraries #######
##########################

library(rgdal)

#####################################
############ Directories ############
#####################################
##create document for butternut pathway
butternut_drive <- "C:\\Users\\eschumacher\\Documents\\GitHub\\butternut"

#####################################
############# Load Files ############
#####################################
##set working directory
setwd(butternut_drive)

##Lon lat files 
butternut_reorg_lonlat <- read.csv("data_files\\after_reorg\\reorg_lon_lat.csv")

##population names document 
butternut_24pop_names <- unique(butternut_reorg_lonlat$Pop)

#############################################
############# Mean Lat and Long #############
#############################################
##calculate mean longitude and latitude for each population
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

##create a minimum and maximum data frame for longitude and latitude 
max_min_df <- matrix(nrow = 2, ncol = 2)
max_min_df[,1] <- c(min_lon, max_lon)
max_min_df[,2] <- c(min_lat, max_lat)
max_min_df<- data.frame(max_min_df)
colnames(max_min_df) <- c("Longitude", "Latitude")
rownames(max_min_df) <- c("Min", "Max")

##write out minimum and maximum coordinates 
write.csv(max_min_df, "data_files\\geographic_files\\max_min_lonlat_df.csv")

###########################
###### Better Pop Map #####
###########################
setwd("data_files\\geographic_files")

##projection you want to use
proj_out <- "+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs"

##WGS84 (lat/long), coordinate format of data downloads
proj_WGS84 <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"

##need North America shapefile
na_shp <- readOGR("bound_p", "boundary_p_v2")
na_shp <- sp::spTransform(na_shp, proj_out)

##north american shapefile 
north_america <- spTransform(na_shp, proj_out)

##Convert to spatial 
coordinates(max_min_df) <- c('Longitude', 'Latitude')
proj4string(max_min_df) <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84")
max_min_df_trans <- spTransform(max_min_df, proj_out)

##output df
max_min_df_trans <- data.frame(max_min_df_trans)

##coordinates of all individuals 
coord_df <- data.frame(butternut_reorg_lonlat$Longitude, butternut_reorg_lonlat$Latitude)
colnames(coord_df) <- c("Longitude", "Latitude")

##convert to spatial data 
coordinates(coord_df) <- c('Longitude', 'Latitude')
proj4string(coord_df) <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84")
coord_df_proj <- spTransform(coord_df, proj_out)

##read in buffer shapefile 
setwd(butternut_drive)
butternut_buffer <- readOGR(dsn = paste0(butternut_drive,"\\data_files\\geographic_files") , 
                            layer = "butternut_buffer")
butternut_buffer_trans <- sp::spTransform(na_shp, proj_out)

##list of 24 colors for populations 
butternut_col <- c("darkmagenta", "darkorchid3", "deeppink3", "darkviolet", "deeppink","lightpink",
                   "lightcoral","orangered3","red4","firebrick1","brown1", "darkgoldenrod1",
                   "darkseagreen", "darkgoldenrod4","gold","deepskyblue", "lightsalmon1","khaki", 
                   "olivedrab4", "olivedrab3","forestgreen","lightblue3","dodgerblue2","orange3")

##add colors to data frame 
butternut_col_df <- matrix(nrow = length(butternut_reorg_lonlat[,1]), ncol = 2)
butternut_col_df[,1] <- as.character(butternut_reorg_lonlat$Pop)
for(i in 1:length(butternut_24pop_names)){
  
  butternut_col_df[butternut_col_df[,1] == paste0(butternut_24pop_names[i]),][,2] <- butternut_col[i]
  
}

##plot mean population coordinates 
setwd("Images\\geographic_images")
pdf("pop_map_allind.pdf", width = 8, height = 8)
plot(north_america, xlim = c(max_min_df_trans[,1]), ylim = c(max_min_df_trans[,2]))
points(coord_df_proj, pch = 17, col = butternut_col_df[,2])
legend('bottomright',legend = butternut_24pop_names, pch = 17, col = butternut_col, cex = 0.6)
dev.off()

##write out coordinate df 
setwd(paste0(butternut_drive, "\\data_files\\geographic_files"))

##combine coordinates and colors by population to create a data frame
butternut_coord_df <- cbind(as.numeric(butternut_mean_lon), as.numeric(butternut_mean_lat), butternut_col)

##name row 
rownames(butternut_coord_df) <- butternut_24pop_names
colnames(butternut_coord_df) <- c("Mean_Lon","Mean_Lat","col")
write.csv(butternut_coord_df, "butternut_coord_df.csv")
