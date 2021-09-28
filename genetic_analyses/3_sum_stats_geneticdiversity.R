#Code below was the first pass at genetic analyses with the increased 
#sampling from butternut's northern population range. 
#These are the next steps in calculating genetic diversity 
#First, we tested linkage disequilibrium, null alleles, 
#and Hardy Weinberg equilibrium.
#Next, a data frame including expected heterozygosity, allelic richness,
#number of alleles, mean longtiude and latitude by population, distance to range edge, and 
#individual numbers. This table is included in full in the supplemental text 
#of this manuscript.

######################################
########### Load Libraries ###########
######################################

library(adegenet)
library(poppr)
library(hierfstat)
library(PopGenReport)
library(pegas)
library(geosphere)
library(sp)
library(rgdal)

##########################################
############# Set directories ############
##########################################
##set directory to all butternut files 
butternut_drive <- "C:\\Users\\eschumacher\\Documents\\GitHub\\Butternut"

##########################################
############ Load in files ###############
##########################################
setwd(butternut_drive)

##convert
#butternutgen_nomd <- arp2gen("Genetic_Analyses\\data_files\\after_reorg\\butternutgen_nomd.arp")

##remove missing data 
butternutgen_reorg <- read.genepop("Genetic_Analyses\\data_files\\after_reorg\\butternutgen_nomd.gen", ncode = 3)

##load lat lon doc for individual names and population names 
butternut_reorg_lonlat <- read.csv("Genetic_Analyses\\data_files\\after_reorg\\reorg_lon_lat.csv")

##load in mean lon and lat document 
butternut_mean_lon_lat <- read.csv("Genetic_Analyses\\data_files\\geographic_files\\butternut_coord_df.csv")

##name individuals in genind doc 
rownames(butternutgen_reorg@tab) <- butternut_reorg_lonlat$Ind

##population names document 
butternut_24pop_names <- unique(butternut_reorg_lonlat$Pop)

##name populations in genind 
levels(butternutgen_reorg@pop) <- butternut_24pop_names

##load range buffer for butternut
butternut_buffer <- readOGR(dsn = paste0(butternut_drive,"\\Genetic_Analyses\\data_files\\geographic_files") , 
                            layer = "butternut_buffer")

############################################################################
####### Run Genetic Diversity Checks like LD, HWE, Null Alleles  ###########
############################################################################

##calculate null alleles 
bn_null_all <- null.all(butternutgen_reorg)

##create null allele table
bn_null_all_df <- matrix(nrow = length(rownames(bn_null_all$null.allele.freq$summary1)),
                         ncol = length(colnames(bn_null_all$null.allele.freq$summary1)))

bn_null_all_df <- bn_null_all$null.allele.freq$summary1

##bn HWE test
bn_hwe <- hw.test(butternutgen_reorg, B = 1000)
bn_hwe_pop <- seppop(butternutgen_reorg) %>% lapply(hw.test, B = 0)

##create table by populations
bn_hwe_bypop <- sapply(bn_hwe_pop, "[", i = TRUE, j = 3)

##name columns
colnames(bn_hwe_bypop) <- butternut_24pop_names

##test for LD
ld_comp <- pair.ia(butternutgen_reorg, sample = 1000)
ld_comp_df <- data.frame(round(ld_comp,digits = 2))

##write out data files 
write.csv(bn_hwe, "Genetic_Analyses\\genetic_analyses_results\\bn_hwe_overall.csv")
write.csv(bn_hwe_bypops, "Genetic_Analyses\\genetic_analyses_results\\bn_hwe_bypop.csv")
write.csv(ld_comp_df, "Genetic_Analyses\\genetic_analyses_results\\ld_loci.csv")

######################################
############# basic stats ############
######################################
##reorg data file 
bn_sumstats <- summary(butternutgen_reorg)

##create poppr file 
BN_poppr <- poppr(butternutgen_reorg)

##expected heterozygosity 
BN_hexp <- BN_poppr[1:24, 10]
##allele numbers by pop 
BN_nall <- bn_sumstats$pop.n.all
##individual numbers
BN_ind <- BN_poppr[1:24, 2:3]
##allelic richness code 
BN_alleles <-bn_sumstats$pop.n.all/length(butternutgen_reorg@loc.n.all)
BN_all_rich <- colMeans(allelic.richness(butternutgen_reorg)$Ar)	

####Geographic analyses for each population
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
##now combine into dataframes 
butternut_coords_df <- cbind(butternut_mean_lon, butternut_mean_lat)[-1,]

#project range extent buffer to the same extent as population lon/lat
butternut_buffer_trans <- sp::spTransform(butternut_buffer, "+proj=longlat +ellps=WGS84 +datum=WGS84")

##now calculate distance to butternut range edge 
butternut_dist <- dist2Line(butternut_coords_df, butternut_buffer_trans)

##add regional names 
butternut_regions <- c("New Brunswick","New Brunswick", "New Brunswick", "New Brunswick", "New Brunswick",
                       "New Brunswick", "Ontario", "Quebec","Ontario","Ontario","Quebec",
                       "United States", "United States", "United States", "United States", "United States",
                       "United States", "United States", "United States", "United States", "United States", 
                       "United States", "United States", "United States")

##create data frame 
butternut_stat_df <- cbind(butternut_coords_df[,c(1:2)], (butternut_dist[,1]/1000), 
                                        BN_ind, BN_nall, BN_all_rich, BN_hexp)

butternut_stat_df <- cbind(butternut_regions, butternut_stat_df)

##name columns and rows 
rownames(butternut_stat_df) <- c(1:24)
colnames(butternut_stat_df) <- c("Region","Mean Longitude", "Mean Latitude", "Distance to Range Edge (km)", 
                                 "Number of Individuals", 
                                 "MLG","Number of Alleles", "Allelic Richness", "Expected Heterozygosity")

##write out csv 
write.csv(butternut_stat_df, "Genetic_Analyses\\genetic_analyses_results\\butternut_stat_df.csv")

##run sessionInfo
sessionInfo()
