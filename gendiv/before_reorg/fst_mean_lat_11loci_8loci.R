
####################
##### Libraries ####
####################

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
library(strataG)
library(plotrix)

###############################################
############## Working Directory ##############
###############################################

setwd("G:/My Drive/Hoban_Lab_Docs/Projects/Butternut_JUCI")

############################################################
########## Genetic Analyses - Mean Lat vs. PWFst ###########
############################################################

##load files 
butternut_24pop_gen <- read.genepop("DataFiles/24Populations/11loci/butternut_24pops.gen", ncode = 3)

##load genotype file
butternut_24pop <- read.csv("DataFiles/24Populations/11loci/butternut_24pop_rel.csv")

##reduce based on missing data 
butternutgen_nomd <- missingno(butternut_24pop_gen, type = "geno", cutoff = 0.25, quiet = FALSE, freq = FALSE)

##data frame to reduce loading
butternut_latlon <- read.csv("DataFiles/24Populations/11loci/butternut_24pops_lonlat.csv")

##first reduce it by missing data 
butternut_latlon_nomd <- butternut_latlon[butternut_latlon$Ind %in% rownames(butternutgen_nomd@tab),]

##calculate PWFst 
butternut_24pop_fst <- as.matrix(pairwise.neifst(butternut_24pop[,-1]))

##code to create name file
butternut_24pop_names <- rownames(butternut_24pop_fst)

##now do lon lat calcs
butternut_mean_lon <- matrix()
butternut_mean_lat <- matrix()

##loops for mean lat/lon
for(pop in butternut_24pop_names){
  
  butternut_mean_lon[pop] <- mean(butternut_latlon[butternut_latlon$Pop == pop,][,3])  
  
}

for(pop in butternut_24pop_names){
  
  butternut_mean_lat[pop] <- mean(butternut_latlon[butternut_latlon$Pop == pop,][,4])  
  
}

##convert to matrix
butternut_mean_lon <- matrix(butternut_mean_lon)
butternut_mean_lat <- matrix(butternut_mean_lat)

##document cleanup
butternut_mean_lon <- butternut_mean_lon[-1]
butternut_mean_lat <- butternut_mean_lat[-1]

#combine into one document for mean long and lat for each pop
butternut_coord_df <- matrix(ncol = 2, nrow = 24)
butternut_coord_df[,1] <- butternut_mean_lon
butternut_coord_df[,2] <- butternut_mean_lat
rownames(butternut_coord_df) <- butternut_24pop_names
colnames(butternut_coord_df) <- c("Mean Lon", "Mean Lat")

#calculate mean PWFst of each population to identify each population's avg difference between all other pops
butternut_meanfst <- matrix()

for(a in 1:24){
  
  butternut_meanfst[a] <- mean(butternut_24pop_fst[a,], na.rm = TRUE)
  
}

##create mean lat and fst table
butternut_fst_lat <- cbind(butternut_meanfst,butternut_mean_lat)
butternut_fst_lat <- data.frame(butternut_fst_lat)

##add names to dataframe
rownames(butternut_fst_lat) <- butternut_24pop_names
colnames(butternut_fst_lat) <- c("Mean PWFst", "Mean Lat")

##Reorganize by geographic range 
butternut_fst_lat <- butternut_fst_lat[c("31","568","1014","7917",
                                         "9101113a","9101113b","151","170","125147","126147",
                                         "171188","APA12PACb","BV", "BVT","BVTGMVT2","CNF", "FKN","GMVT1", "MC",    
                                         "OZ12","SFNF12", "WHWI", "WWI","VSH123"),]
##Create color column
butternut_fst_lat$Col <- NA

butternut_fst_lat[1:6,3] <- "firebrick1"
butternut_fst_lat[7:11,3] <- "firebrick4"
butternut_fst_lat[12:24,3] <- "dodgerblue"

##create linear model for plotting comparison with mean lat compared to mean population fst
butternut_24pop_latfst_lm <- lm(as.numeric(butternut_fst_lat[,1])~as.numeric(butternut_fst_lat[,2]))
butternut_24pop_latfst_lm_sum <- summary(butternut_24pop_latfst_lm)

##create r squared and p-value legend
r2_24pop_latfst <- butternut_24pop_latfst_lm_sum$adj.r.squared
coef_24pop_latfst <- butternut_24pop_latfst_lm_sum$coefficients
pvalue_24pop_latfst <- coef_24pop_latfst[2,4]
rp_24pop_latfst <- vector('expression',2)
rp_24pop_latfst[1] <- substitute(expression(italic(R)^2 == MYVALUE), 
                                 list(MYVALUE = format(r2_24pop_latfst,dig=3)))[2]
rp_24pop_latfst[2] <- substitute(expression(italic(p) == MYOTHERVALUE), 
                                 list(MYOTHERVALUE = format(pvalue_24pop_latfst, digits = 2)))[2]

##butternut poppr document
butternut_24pop_poppr <- poppr(butternutgen_nomd)
butternut_24pop_ind <- butternut_24pop_poppr[1:24,2]

##plot mean lat vs. mean pwfst 
pdf(file="Graphical_Stat_Results/PostIndRemoval/24pop/11loci/latfst_24pops_11loci.pdf", width=10,height=6)
plot(as.numeric(butternut_fst_lat[,1])~as.numeric(butternut_fst_lat[,2]), pch = 16,  cex = c(butternut_24pop_ind/50), xlab = "Mean Latitude", ylab = "Mean Fst", col = butternut_fst_lat[,3], ylim = c(0,0.04))
text(as.numeric(butternut_fst_lat[,1])~as.numeric(butternut_fst_lat[,2]), label = butternut_24pop_names, cex = 0.8, pos = 3,offset = 0.5)
legend('topleft', legend = rp_24pop_latfst, bty = 'n', border = "black", pt.cex = 1, cex = 0.8)
legend('top', legend = c("NB", "Ontario","US"), pch = 16, col = c("firebrick1","firebrick4","dodgerblue"))
abline(butternut_24pop_latfst_lm, col = "darkblue")
dev.off()

###############################################
########## Loci Reduction - 8 Loci ############
###############################################

##remove loci - B157, B212_1, B262
butternutgen_24pop_8loci <- butternutgen_nomd[loc=-c(5,6,10)]

###############################################
######## 8 Loci - Genetic Diversity ###########
###############################################
###set up long/lat document
butternut_24pop_8loci_lonlat <- butternut_latlon_nomd[,-c(butternut_latlon_nomd$B157_1,butternut_latlon_nomd$B157_2,
                                                         butternut_latlon_nomd$B212_1, butternut_latlon_nomd$B212_2,
                                                         butternut_latlon_nomd$B262_1, butternut_latlon_nomd$B262_2)]

##Lat vs. PWFst
butternut_24pop_8loci_fst <- as.matrix(pairwise.neifst(butternut_24pop_8loci_lonlat[,-1]))

##calculating mean fst
butternut_8loci_meanfst <- matrix()

for(a in 1:24){
  
  butternut_8loci_meanfst[a] <- mean(butternut_24pop_8loci_fst[a,], na.rm = TRUE)
  
}

butternut_8loci_fst_lat <- cbind(butternut_coord_df[,2],butternut_8loci_meanfst)
rownames(butternut_8loci_fst_lat) <- butternut_24pop_names
colnames(butternut_8loci_fst_lat) <- c("Mean Lat","Mean Fst")

##reorganize the document 
fst_lat_8loci <- butternut_8loci_fst_lat[c("31","568","1014","7917",
                                         "9101113a","9101113b","151","170","125147","126147",
                                         "171188","APA12PACb","BV", "BVT","BVTGMVT2","CNF", "FKN","GMVT1", "MC",    
                                         "OZ12","SFNF12", "WHWI", "WWI","VSH123"),]

##create a dataframe 
fst_lat_8loci <- data.frame(fst_lat_8loci)

##add colors 
fst_lat_8loci$Col <- butternut_fst_lat$Col

##now plot the comparison
fst_lat_8loci_lm <- lm(fst_lat_8loci[,2]~fst_lat_8loci[,1])
fst_lat_8loci_summary <- summary(fst_lat_8loci_lm)

r2_8loci_fst <- fst_lat_8loci_summary$adj.r.squared
coef_8loci_fst <- fst_lat_8loci_summary$coefficients
pvalue_8loci_fst <- coef_8loci_fst[2,4]
rp_8loci_fst <- vector('expression',2)
rp_8loci_fst[1] <- substitute(expression(italic(R)^2 == MYVALUE), 
                              list(MYVALUE = format(r2_8loci_fst,dig=3)))[2]
rp_8loci_fst[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                             list(MYOTHERVALUE = format(pvalue_8loci_fst, digits = 2)))[2]

pdf(file="Graphical_Stat_Results/PostIndRemoval/24pop/8loci/fst_lat_8loci.pdf", width=6,height=6)
plot(fst_lat_8loci[,1],fst_lat_8loci[,2], xlab = "Mean Lat", ylab = "Mean Population Fst", pch = 16, cex = c(butternut_24pop_ind/50), col = fst_lat_8loci[,3])
legend('topleft', legend = rp_8loci_fst, bty = 'n', border = "black", pt.cex = 1, cex = 0.8)
text(butternut_8loci_fst_lat[,1],butternut_8loci_fst_lat[,2], label = rownames(butternut_8loci_fst_lat), pos = 2.5, offset = 0.5, cex = 0.8)
abline(fst_lat_8loci_lm, col = "blue")
legend('top', legend = c("NB", "Ontario","US"), pch = 16, col = c("firebrick1","firebrick4","dodgerblue"))
dev.off()

##write out data files 
write.csv(butternut_fst_lat, "DataFiles\\24Populations\\11loci\\fst_lat_11loci.csv")
write.csv(butternut_24pop_8loci_lonlat, "DataFiles\\24Populations\\8loci\\lon_lat_8loci.csv")
write.csv(fst_lat_8loci, "DataFiles\\24Populations\\8loci\\fst_lat_8loci.csv")
write.csv(butternut_coord_df,"DataFiles\\24Populations\\mean_coords.csv")

