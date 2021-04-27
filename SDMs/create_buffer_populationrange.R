##########################
######## Libraries #######
##########################



#####################################
############ Directories ############
#####################################
butternut_drive <- "G:\\My Drive\\Hoban_Lab_Docs\\Projects\\Butternut_JUCI"

##set wd
setwd(butternut_drive)

######load in presence records
butternut_range <- read.csv("SDMs\\worldclim_only\\InputFiles\\occurrence_noauto_noproj.csv")
butternut_range2 <- butternut_range[,-1]

#####calculate maximums and minimums
min_lon <- min(butternut_range$Longitude)
max_lon <- max(butternut_range$Longitude)
min_lat <- min(butternut_range$Latitude)
max_lat <- max(butternut_range$Latitude)

######make spatial points
coordinates(butternut_range2) <- c('Longitude', 'Latitude')
proj4string(butternut_range2) <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84")

#########################################
########## Calculate Buffer #############
#########################################
setwd(butternut_drive)
######load in presence records
butternut_occ <- read.csv("SDMs\\worldclim_only\\InputFiles\\occurrence_noauto_noproj.csv")
butternut_occ <- butternut_range[,-1]

#projection you want to use
proj_out <- "+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs"

# WGS84 (lat/long), coordinate format of data downloads
proj_WGS84 <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"

##set working direcory
setwd("Graphical_Stat_Results\\PostIndRemoval\\GeographicImages")
##load in North America shapefile
na_shp <- readOGR("bound_p", "boundary_p_v2")
na_shp <- sp::spTransform(na_shp, proj_out)

# make spatial points dataframe using occurrence data
coords <- butternut_occ %>% dplyr::select("Longitude","Latitude")
names(coords) <- c("x","y")
coordinates(coords) <- coords
sp::proj4string(coords) <- proj_WGS84
coords_t <- sp::spTransform(coords, proj_out, multiline = "NO")
coords <- coordinates(coords_t)
spdf <- SpatialPointsDataFrame(coords = coords, data = butternut_occ, proj4string = CRS(proj_out))

# crop occurrence points to area of interest (ie, North America)
crop <- crop(x = spdf, y = na_shp@bbox)

# create 100km buffer around occurrence points
butternut_buff <- gBuffer(crop, width = 100000)

##write out buffer file 
shapefile(butternut_buff,"butternut_buffer")
