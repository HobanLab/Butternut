######################################
########### Load Libraries ###########
######################################

library(adegenet)
library(poppr)
library(pegas)
library(PopGenReport)
library(hierfstat)

#####################################
############# Load Files ############
#####################################
##set working directory
butternut_dir <- "G:/My Drive/Hoban_Lab_Docs/Projects/Butternut_JUCI"
setwd(butternut_dir)

##load reorg genind
butternutgen_reorg <- read.genepop("DataFiles/24Populations/reorg/reorg_gen_24pop.gen", ncode = 3)

##load lat lon doc 
butternut_reorg_lonlat <- read.csv("DataFiles/24Populations/reorg/reorg_lon_lat.csv")

##create name doc
butternut_24pop_names <- unique(butternut_reorg_lonlat$Pop)

##name genind doc 
rownames(butternutgen_reorg@tab) <- butternut_reorg_lonlat$Ind

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
BN_alleles<-bn_sumstats$pop.n.all/length(butternutgen_reorg@loc.n.all)
BN_all_rich<-colSums(allelic.richness(butternutgen_reorg)$Ar)/length(butternutgen_reorg@loc.n.all)	#remember to divide by number of loci!!

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

##create pop names 
butternut_24pop_names <- unique(butternut_reorg_lonlat$Pop)

##create data frame 
butternut_stat_df <- signif(cbind(butternut_mean_lon, butternut_mean_lat, BN_ind, BN_nall, BN_all_rich, BN_hexp),3)

##name columns and rows 
rownames(butternut_stat_df) <- butternut_24pop_names
colnames(butternut_stat_df) <- c("Mean Longitude", "Mean Latitude", "Number of Individuals", "MLG","Number of Alleles", "Allelic Richness", "Expected Heterozygosity")

##write out csv 
write.csv(butternut_stat_df, "Graphical_Stat_Results\\PostIndRemoval\\24pop\\Reorg_Results\\GeneticDiversity\\butternut_stat_df.csv")

##############################
####### Run Checks ###########
##############################

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
bn_hwe_reorg <- sapply(bn_hwe_pop, "[", i = TRUE, j = 3)

##name columns
colnames(bn_hwe_reorg) <- butternut_24pop_names

##test for LD
ld_comp <- pair.ia(butternutgen_reorg, sample = 1000)
ld_comp_df <- data.frame(round(ld_comp,digits = 2))

##write out data files 
write.csv(bn_hwe, "Graphical_Stat_Results\\PostIndRemoval\\24pop\\Reorg_Results\\GeneticDiversity\\bn_hwe.csv")
write.csv(bn_hwe_reorg, "Graphical_Stat_Results\\PostIndRemoval\\24pop\\Reorg_Results\\GeneticDiversity\\bn_hwe_reorg.csv")
write.csv(ld_comp_df, "Graphical_Stat_Results\\PostIndRemoval\\24pop\\Reorg_Results\\GeneticDiversity\\ld_loci.csv")

