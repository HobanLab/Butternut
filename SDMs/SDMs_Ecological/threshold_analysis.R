
library(dplyr)



##Load in presence and absence points 
butternut_pres <- butternut_sel_df[butternut_sel_df$PA == 1,]
butternut_abs <- butternut_sel_df[butternut_sel_df$PA == 0,]

##create points 
coordinates(butternut_pres) <- c('Longitude', 'Latitude')
proj4string(butternut_pres) <- CRS(projection)

coordinates(butternut_abs) <- c('Longitude', 'Latitude')
proj4string(butternut_abs) <- CRS(projection)


##extract to points 
hsm_p <- extract(butternut_prediction_map, butternut_pres)
hsm_a <- extract(butternut_prediction_map, butternut_abs)

##extract df
hsm_p_df <- data.frame(hsm_p)
hsm_a_df <- data.frame(hsm_a)

##cbind coords 
hsm_p_df <- data.frame(butternut_sel_df[butternut_sel_df$PA == 1,], hsm_p)
hsm_p_df <- hsm_p_df[,-c(3:8)]

##calculate 95% of points
hsm_p_remove <- round(length(hsm_p_df$hsm_p)*0.05 )

##
new_hsm_p <- hsm_p_df %>% 
  # desc orders from largest to smallest
  arrange(desc(hsm_p))

##calculate which rows to remove 

hsm_r_upper <- length(hsm_p_df$hsm_p) - hsm_p_remove 

##
suitability_hsm_df <- new_hsm_p[-c(hsm_r_upper:length(new_hsm_p$hsm_p)),]
##last row

##lowest end of suitability = 0.2425463
hsm_scaled <- butternut_prediction_map*1.9
par(mfrow=c(1,1))
hsm_prob = (exp(hsm_scaled) / (1+ exp(hsm_scaled)))
hsm_prob=(exp(butternut_prediction_map)/(1+exp(butternut_prediction_map)))
plot(hsm_prob)
points(butternut_pres)

new_ext <- extract(hsm_prob, butternut_pres)

new_df_ext <- data.frame(new_ext)

##
e <- evaluate(butternut_pres, butternut_abs, butternut_model, butternut_model_stack)
tr <- threshold(e, 'spec_sens')


plot(butternut_prediction_map > tr)

###with new threshold 
setwd("G:\\Shared drives\\Emily_Schumacher\\SDMs\\worldclim_only_SDM\\HSM_rasters_dif_timepoints")

###load rasters 
hsm_list <- list.files(pattern = "_hsm.tif$")
hsm_raster_list <- list()

for(r in 1:length(hsm_list)){
  
  hsm_raster_list[[r]] <- raster(hsm_list[[r]])
  
  
}

##set up time period list
time_periods <- c("ba - 14.7-12.9 ka", 
                  "eh - 11.7-8.326 ka", 
                  "hs - 17.0-14.7 ka", 
                  "lgm - 22 ka",
                  "lh - 4.2-0.3 ka",
                  "lig - 130 ka",
                  "mid_hol - 6 ka",
                  "present - 1970 - 2000", 
                  "yds - 12.9-11.7 ka")


##
setwd(paste0(worldclim_only, "tpr_tnr_maps"))

pdf("ba_tr.pdf", width = 8, height = 6)
plot(hsm_raster_list[[1]] > tr, main = paste0(time_periods[[1]]), legend = FALSE)
dev.off()

pdf("eh_tr.pdf", width = 8, height = 6)
plot(hsm_raster_list[[2]] > tr, main = paste0(time_periods[[2]]), legend = FALSE)
dev.off()

pdf("hs_tr.pdf", width = 8, height = 6)
plot(hsm_raster_list[[3]] > tr, main = paste0(time_periods[[3]]), legend = FALSE)
dev.off()

pdf("lgm_tr.pdf", width = 8, height = 6)
plot(hsm_raster_list[[4]] > tr, main = paste0(time_periods[[4]]), legend = FALSE)
dev.off()

pdf("lh_tr.pdf", width = 8, height = 6)
plot(hsm_raster_list[[5]] > tr, main = paste0(time_periods[[5]]), legend = FALSE)
dev.off()

pdf("lig_tr.pdf", width = 8, height = 6)
plot(hsm_raster_list[[6]] > tr, main = paste0(time_periods[[6]]), legend = FALSE)
dev.off()

pdf("mid_hol_tr.pdf", width = 8, height = 6)
plot(hsm_raster_list[[7]] > tr, main = paste0(time_periods[[7]]), legend = FALSE)
dev.off()

pdf("present_tr.pdf", width = 8, height = 6)
plot(hsm_raster_list[[8]] > tr, main = paste0(time_periods[[8]]), legend = FALSE)
dev.off()

pdf("yds_tr.pdf", width = 8, height = 6)
plot(hsm_raster_list[[9]] > tr, main = paste0(time_periods[[9]]), legend = FALSE)
dev.off()



