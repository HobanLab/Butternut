
######################################################################
########################### Libraries ################################
######################################################################

library(raster)
library(sp)
library(sf)
library(rworldmap)
library(rgdal)
library(spdep)
library(rgeos)
library(dismo)
library(gbm)
library(AUC)
library(ggplot2)
library(plyr)
library(HH)
library(ltm)
library(geosphere)

#########################################################################
##################### Load Files ########################################
#########################################################################
setwd("G:\\Shared drives\\Emily_Schumacher\\SDMs\\worldclim_only_SDM")

projection <- c("+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs")


###load rasters 
present_hsm <- raster("present_hsm_project.tif")
lgm_hsm <- raster("lgm_hsm_project.tif")
mid_hol_hsm <- raster("mid_hol_hsm_project.tif")
lig_hsm <- raster("lig_hsm_project.tif")

###project rasters 

present_hsm_proj <- projectRaster(present_hsm, crs = projection)
lgm_hsm_proj <- projectRaster(lgm_hsm, crs = projection)
mid_hol_hsm_proj <- projectRaster(mid_hol_hsm, crs = projection)
lig_hsm_proj <- projectRaster(lig_hsm, crs = projection)

###convert to presence rasters 
present_pres_rast <- present_hsm_proj > 0.63
lgm_pres_rast <- lgm_hsm_proj > 0.63
mid_hol_rast <- mid_hol_hsm_proj > 0.63
lig_pres_rast <- lig_hsm_proj > 0.63

##write out as a plot

##setwd
setwd("G:\\Shared drives\\Emily_Schumacher\\butternut_publication_figures\\dist_edge")

##need North America shapefile
na_shp <- readOGR("bound_p", "boundary_p_v2")

na_shp_proj <- spTransform(na_shp, CRS(projection))

###plot this 
setwd("G:\\Shared drives\\Emily_Schumacher\\SDMs\\worldclim_only_SDM")

pdf("presentday_hsm.pdf")
plot(present_pres_rast)
plot(na_shp_proj, add = TRUE)
dev.off()

##only presence
pres_only_rast <- present_pres_rast$layer > 0

##convert to shp
pres_poly <- rasterToPolygons(present_pres_rast)
lgm_poly <- rasterToPolygons(lgm_pres_rast)
mid_hol_poly <- rasterToPolygons(mid_hol_rast)
lig_poly <- rasterToPolygons(lig_pres_rast)

##now create only presence areas 
pres_poly_only <- pres_poly[pres_poly$layer == 1,]
lgm_pres_only <- lgm_poly[lgm_poly$layer == 1,]
mid_hol_pres_only <- mid_hol_poly[mid_hol_poly$layer == 1,]
lig_pres_only <- lig_poly[lig_poly$layer == 1,]

##now project
pres_poly_only_proj <- spTransform(pres_poly_only, 
                                   CRS(projection))

lgm_pres_only_proj <- spTransform(lgm_pres_only, 
                                  CRS(projection))

mid_hol_pres_only_proj <- spTransform(mid_hol_pres_only, 
                                      CRS(projection))

lig_pres_only_proj <- spTransform(lig_pres_only, 
                                  CRS(projection))

##calculate area
pres_area <- gArea(pres_poly_only_proj)/1000
lgm_area <- gArea(lgm_pres_only_proj)/1000
mid_hol_area <- gArea(mid_hol_pres_only_proj)/1000
lig_area <- gArea(lig_pres_only_proj)/1000

##write table 
area_df <- matrix(nrow = 4, ncol = 2)
area_df[,1] <- c("LIG", "LGM", "Mid Hol", "Present")
area_df[,2] <- c(lig_area, lgm_area, mid_hol_area, pres_area)

write.csv(area_df, "suitable_habitat_area_df.csv")

##intersect different areas 
pres_lgm_intersect <- gIntersection(pres_poly_only_proj, 
                                    lgm_pres_only_proj, 
                                    byid = FALSE)

lgm_lig_intersect <- gIntersection(lgm_pres_only_proj, 
                                   lig_pres_only_proj, 
                                   byid = FALSE)

lgm_mid_hol_intersect <- gIntersection(lgm_pres_only_proj, 
                                       mid_hol_pres_only_proj, 
                                       byid = FALSE)

lig_mid_hol_intersect <- gIntersection(lig_pres_only_proj, 
                                       mid_hol_pres_only_proj, 
                                       byid = FALSE)

pres_mid_hol_intersect <- gIntersection(pres_poly_only_proj, 
                                        mid_hol_pres_only_proj, 
                                        byid = FALSE)

pres_lig_intersect <- gIntersection(pres_poly_only_proj, 
                                    lig_pres_only_proj, 
                                    byid = FALSE)

###calculate area at each point 

pres_lgm_intersect_area <- gArea(pres_lgm_intersect)/1000

lgm_lig_intersect_area <- gArea(lgm_lig_intersect)/1000
+
lgm_mid_hol_intersect_area <- gArea(lgm_mid_hol_intersect)/1000

lig_mid_hol_intersect_area <- gArea(lig_mid_hol_intersect)/1000

pres_lig_intersect_area <- gArea(pres_lig_intersect)/1000

pres_mid_hol_intersect_area <- gArea(pres_mid_hol_intersect)/1000
