##########################
######## Libraries #######
##########################

library(adegenet)
library(hierfstat)
library(poppr)
library(PopGenReport)

#####################################
############ Directories ############
#####################################
butternut_drive <- "G:\\My Drive\\Hoban_Lab_Docs\\Projects\\Butternut_JUCI"

#####################################
############# Load Files ############
#####################################
setwd(butternut_drive)

##load current working reorganized genepop file 
butternutgen_reorg <- read.genepop("DataFiles\\24Populations\\reorg\\reorg_gen_24pop.gen", ncode = 3)

##LOAD GENALEX CSV
reorg_genalex <- read.csv("DataFiles\\24Populations\\reorg\\reorg_genalex.csv")
ind_names <- unique(reorg_genalex[3:1637,1])

rownames(butternutgen_reorg@tab) <- ind_names

##GESTE Results
reorg_geste_fst <- read.csv("DataFiles\\24Populations\\reorg\\GESTE_fst.csv")

##Lon lat files 
butternut_reorg_lonlat <- read.csv("DataFiles\\24Populations\\reorg\\reorg_lon_lat.csv")

##create population name file 
butternut_24pop_names <- unique(butternut_reorg_lonlat$Pop)

##calculate poppr
butternut_24pop_poppr <- poppr(butternutgen_reorg)

#####################################
######### Comparison Plot ###########
#####################################
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

##fst and mean latitude
fst_lat <- data.frame(butternut_mean_lat,reorg_geste_fst[,2:5])
rownames(fst_lat) <- butternut_24pop_names
fst_lat$Col <- NA
fst_lat[1:6,6] <- "firebrick1"
fst_lat[7:11,6] <- "firebrick4"
fst_lat[12:24,6] <- "dodgerblue"

######regression
fst_lat_lm <- lm(fst_lat[,2]~fst_lat[,1])
fst_lat_lm_sum <- summary(fst_lat_lm)

#####create r2 and p-value
##create r squared and p-value legend
latfst_r2 <- fst_lat_lm_sum$adj.r.squared
latfst_coef <- fst_lat_lm_sum$coefficients
latfst_pvalue <- latfst_coef[2,4]
latfst_rp <- vector('expression',2)
latfst_rp[1] <- substitute(expression(italic(R)^2 == MYVALUE), 
                                 list(MYVALUE = format(latfst_r2,dig=3)))[2]
latfst_rp[2] <- substitute(expression(italic(p) == MYOTHERVALUE), 
                                 list(MYOTHERVALUE = format(latfst_pvalue, digits = 2)))[2]

##plot against each other
pdf("G:\\My Drive\\Hoban_Lab_Docs\\Projects\\Butternut_JUCI\\Graphical_Stat_Results\\PostIndRemoval\\24pop\\Reorg_Results\\STR\\geste_fst_linear.pdf", width = 10, height = 8)
##plot fst lat in comparison
plot(fst_lat$Mean~fst_lat$butternut_mean_lat, ylim = c(0,0.15), pch = 17, xlab = "Mean Latitude", 
       ylab = "GESTE Fst", col = fst_lat[,6], cex = (butternut_24pop_poppr[1:24,2]/50), 
       main = "GESTE Fst Comapred to Mean Latitude")
##add population names
text(fst_lat[,1], fst_lat[,2], labels = butternut_24pop_names, cex = 0.8, pos = 3)
##add model
abline(fst_lat_lm, col = "dodgerblue4")
legend('topleft', legend = latfst_rp, bty = 'n', border = "black", pt.cex = 1, cex = 0.8)
legend('top', legend = c("New Brunswick", "Ontario", "United States"), pch = 17, col = c("firebrick1", "firebrick4","dodgerblue"))
dev.off()
