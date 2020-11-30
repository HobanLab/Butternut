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
setwd("G:\\Shared drives\\Emily_Schumacher\\SDMs\\worldclim_only_SDM\\HSM_rasters_dif_timepoints")

projection <- c("+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs")

###load rasters 
hsm_list <- list.files(pattern = "_hsm.tif$")
hsm_raster_list <- list()
hsm_suitable_list <- list()

for(r in 1:length(hsm_list)){
  
  hsm_raster_list[[r]] <- raster(hsm_list[[r]])
  
  hsm_suitable_list[[r]] <- hsm_raster_list[[r]] > 0.63
  
}

##convert to all to shapefiles 
hsm_poly <- list()

for(p in 1:length(hsm_suitable_list)){
  
  hsm_poly[[p]] <- rasterToPolygons(hsm_suitable_list[[p]])
  
}


##now create only presence areas 
pres_only_poly <- list()

for(pres in 1:length(hsm_poly)){
  
  pres_only_poly[[pres]] <- hsm_poly[[pres]][hsm_poly[[pres]]$layer == 1,]
  
}


##calculate area
area_list <- list()
area_df <- matrix(nrow = length(pres_only_poly), ncol = 2)

for(a in 1:length(pres_only_poly)){
  
  area_list[[a]] <- gArea(pres_only_poly[[a]])/1000
  
  area_df[a,1] <- names(hsm_raster_list[[a]])
  
  area_df[a,2] <- area_list[[a]]
  
}
setwd("G:\\Shared drives\\Emily_Schumacher\\SDMs\\worldclim_only_SDM")

##write table 
area_df <- data.frame(area_df)
area_df$time_frame <- c("14700_12900", "11700_8326","17000_14700",
                        "22000", "4200_300","130000","6000","0",
                        "12900_11700")
##set up data frame
colnames(area_df) <- c("Time Period Name","Area of Suitable Habitat",
                       "YBP")

write.csv(area_df, "suitable_habitat_area.csv")

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


##intersection matrix? 
ba_intersect_list <- list()
eh_intersect_list <- list()
hs_intersect_list <- list()
lgm_intersect_list <- list()
lh_intersect_list <- list()
lig_intersect_list <- list()
mid_hol_intersect_list <- list()
present_intersect_list <- list()
yds_intersect_list <- list()


for(i in 1:length(pres_only_poly)){

  ba_intersect_list[[i]] <- gIntersection(pres_only_poly[[1]], pres_only_poly[[i]])
  
}

for(i in 1:length(pres_only_poly)){
  
  eh_intersect_list[[i]] <- gIntersection(pres_only_poly[[2]], pres_only_poly[[i]])
  
}

for(i in 1:length(pres_only_poly)){
  hs_intersect_list[[i]] <- gIntersection(pres_only_poly[[3]], pres_only_poly[[i]])
  
}

for(i in 1:length(pres_only_poly)){
  lgm_intersect_list[[i]] <- gIntersection(pres_only_poly[[4]], pres_only_poly[[i]])
  
}
for(i in 1:length(pres_only_poly)){
  
  lh_intersect_list[[i]] <- gIntersection(pres_only_poly[[5]], pres_only_poly[[i]])
  
}

for(i in 1:length(pres_only_poly)){
  lig_intersect_list[[i]] <- gIntersection(pres_only_poly[[6]], pres_only_poly[[i]])
  
}

for(i in 1:length(pres_only_poly)){
  mid_hol_intersect_list[[i]] <- gIntersection(pres_only_poly[[7]], pres_only_poly[[i]])
  
}

for(i in 1:length(pres_only_poly)){
  present_intersect_list[[i]] <- gIntersection(pres_only_poly[[8]], pres_only_poly[[i]])
  
}

for(i in 1:length(pres_only_poly)){
  yds_intersect_list[[i]] <- gIntersection(pres_only_poly[[9]], pres_only_poly[[i]])
  
}


##calculate area at each point 
intersect_df <- matrix(nrow = 9, ncol = 9)
time_periods <- c("ba", "eh", "hs", "lgm","lh","lig","mid_hol",
                  "present", "yds")
rownames(intersect_df) <- time_periods
colnames(intersect_df) <- time_periods

##now have my loops 
for(a in 1:length(ba_intersect_list)){
  intersect_df[a,1] <- gArea(ba_intersect_list[[a]])/1000
  
}

for(a in 1:length(eh_intersect_list)){
  intersect_df[a,2] <- gArea(eh_intersect_list[[a]])/1000
  
}

for(a in 1:length(hs_intersect_list)){
  intersect_df[a,3] <- gArea(hs_intersect_list[[a]])/1000
  
}

for(a in 1:length(lgm_intersect_list)){
  intersect_df[a,4] <- gArea(lgm_intersect_list[[a]])/1000
  
}

for(a in 1:length(lh_intersect_list)){
  intersect_df[a,5] <- gArea(lh_intersect_list[[a]])/1000
  
}

for(a in 1:length(lig_intersect_list)){
  intersect_df[a,6] <- gArea(lig_intersect_list[[a]])/1000
  
}

for(a in 1:length(mid_hol_intersect_list)){
  intersect_df[a,7]<- gArea(mid_hol_intersect_list[[a]])/1000
  
}

for(a in 1:length(present_intersect_list)){
  intersect_df[a,8] <- gArea(present_intersect_list[[a]])/1000
  
}

for(a in 1:length(yds_intersect_list)){
  intersect_df[a,9] <- gArea(yds_intersect_list[[a]])/1000
  
}

##reorder columns
intersect_df <- intersect_df[,c("lig", "lgm","hs","ba","yds","eh",
                               "mid_hol","lh","present")]
intersect_df <- intersect_df[c("lig", "lgm","hs","ba","yds","eh",
                                "mid_hol","lh","present"),]

###calculate % overlap of intersection 

intersect_percent <- matrix(nrow = 9, ncol = 9)

for(i in 1:9){
  for(a in 1:9){
    intersect_percent[a,i] <- (intersect_df[a,i]/intersect_df[a,a])*100
    
  }
  
}

rownames(intersect_percent) <- time_periods
colnames(intersect_percent) <- time_periods

write.csv(intersect_percent, "intersect_percent.csv")




