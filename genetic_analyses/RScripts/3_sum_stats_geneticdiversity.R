#Code below was the first pass at genetic analyses with the increased 
#sampling from butternut's northern population range. 
#These are the next steps in calculating genetic diversity 
#First, we tested linkage disequilibrium, null alleles, 
#and Hardy Weinberg equilibrium.
#Next, a data frame including expected heterozygosity, allelic richness,
#number of alleles, mean longtiude and latitude by population, distance to range edge, and 
#individual numbers. This table is included in full in the supplemental text 
#of this manuscript.

###########################
#     Load Libraries      #
###########################

library(adegenet)
library(poppr)
library(hierfstat)
library(PopGenReport)
library(pegas)
library(geosphere)
library(sp)
library(rgdal)

#########################
#     Load in files     #
#########################
setwd("../data_files/after_reorg")

##comment in if you to convert from .arp to .gen file 
#butternutgen_nomd <- arp2gen("Genetic_Analyses\\data_files\\after_reorg\\butternutgen_nomd.arp")

#load 24 pop genepop file as a genind object 
butternutgen_reorg <- read.genepop("butternut_24pop_nomd.gen", ncode = 3)

#load lat lon doc for individual names and population names 
butternut_reorg_lonlat <- read.csv("butternut_24pop_lonlat.csv")

##load in mean lon and lat document 
butternut_mean_lonlat <- read.csv("../geographic_files/butternut_coord_df.csv")

#name individuals in genind doc 
rownames(butternutgen_reorg@tab) <- butternut_reorg_lonlat$Ind

#create pop name object
butternut_24pop_names <- unique(butternut_reorg_lonlat$Pop)

#name populations in genind 
levels(butternutgen_reorg@pop) <- butternut_24pop_names

#load range buffer for butternut
butternut_buffer <- readOGR(dsn = "../geographic_files" , 
                            layer = "butternut_buffer")

###########################################################################
#     Check for Null Alleles, HWE Deviations, Linkage Disequilibrium      #
###########################################################################
#calculate null alleles 
butternut_nullall <- null.all(butternutgen_reorg)

#create null allele table
butternut_nullall_df <- matrix(nrow = length(rownames(butternut_nullall$null.allele.freq$summary1)),
                         ncol = length(colnames(butternut_nullall$null.allele.freq$summary1)))

#save results in data frame
butternut_nullall_df <- bn_null_all$null.allele.freq$summary1

#bn HWE test
butternut_hwe_pop <- seppop(butternutgen_reorg) %>% lapply(hw.test, B = 0)

#create table by populations
butternut_hwe_pop_df <- sapply(butternut_hwe_pop, "[", i = TRUE, j = 3)

##name columns
colnames(butternut_hwe_pop_df) <- butternut_24pop_names

#test forLinkage Disequilibrium 
butternut_ld <- pair.ia(butternutgen_reorg, sample = 1000)
butternut_ld_df <- data.frame(round(butternut_ld,digits = 2))

#write out data files 
write.csv(bn_hwe, "../../genetic_analyses_results/Diversity_Analyses/bn_hwe_overall.csv")
write.csv(bn_hwe_bypops, "Genetic_Analyses\\genetic_analyses_results\\bn_hwe_bypop.csv")
write.csv(ld_comp_df, "Genetic_Analyses\\genetic_analyses_results\\ld_loci.csv")

###################################
#     Genetic Diversity Stats     #
###################################
#reorg data file 
butternut_sumstats <- summary(butternutgen_reorg)

#create poppr file 
butternut_poppr <- poppr(butternutgen_reorg)

#expected heterozygosity 
butternut_hexp <- butternut_poppr[1:24, 10]

#allele numbers by pop 
butternut_nall <- butternut_sumstats$pop.n.all

#individual numbers
butternut_ind <- butternut_poppr[1:24, 2:3]

#allelic richness code 
butternut_all <- butternut_sumstats$pop.n.all/length(butternutgen_reorg@loc.n.all)
butternut_allrich <- colMeans(allelic.richness(butternutgen_reorg)$Ar)	

##Geographic analyses for each population
#calculate mean longitude and latitude for each population
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

#add regional names 
butternut_regions <- c("New Brunswick","New Brunswick", "New Brunswick", "New Brunswick", "New Brunswick",
                       "New Brunswick", "Ontario", "Quebec","Ontario","Ontario","Quebec",
                       "United States", "United States", "United States", "United States", "United States",
                       "United States", "United States", "United States", "United States", "United States", 
                       "United States", "United States", "United States")

#create summary stats df
butternut_gendiv_stat_df <- cbind(butternut_coords_df[,c(1:2)], (butternut_dist[,1]/1000), 
                           butternut_ind, butternut_nall, butternut_allrich, butternut_hexp)

butternut_gendiv_stat_df <- cbind(butternut_regions, butternut_gendiv_stat_df)

#name columns and rows 
rownames(butternut_gendiv_stat_df) <- c(1:24)
colnames(butternut_gendiv_stat_df) <- c("Region","Mean Longitude", "Mean Latitude", "Distance to Range Edge (km)", 
                                 "Number of Individuals", 
                                 "MLG","Number of Alleles", "Allelic Richness", "Expected Heterozygosity")

#write out csv 
write.csv(butternut_gendiv_stat_df, "../../genetic_analyses_results/Diversity_Analyses/butternut_gendiv_stat_df.csv")

##run sessionInfo
sessionInfo()
