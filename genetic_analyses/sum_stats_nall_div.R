######################################
########### Load Libraries ###########
######################################

library(adegenet)
library(poppr)
library(hierfstat)

##########################################
############# Set directories ############
##########################################
##set directory to all butternut files 
butternut_drive <- "C:\\Users\\eschumacher\\Documents\\GitHub\\butternut"

##########################################
############ Load in files ###############
##########################################
setwd(butternut_drive)

##load reorg genind
butternutgen_reorg <- read.genepop("data_files\\after_reorg\\reorg_gen_24pop.gen", ncode = 3)

##load lat lon doc for individual names and population names 
butternut_reorg_lonlat <- read.csv("data_files\\after_reorg\\reorg_lon_lat.csv")

##load in mean lon and lat document 
butternut_mean_lon_lat <- read.csv("data_files\\geographic_files\\butternut_coord_df.csv")

##create population name doc
butternut_24pop_names <- unique(butternut_reorg_lonlat$Pop)

##name individuals in genind doc 
rownames(butternutgen_reorg@tab) <- butternut_reorg_lonlat$Ind
##name populations in genind doc 
levels(butternutgen_reorg@pop) <- butternut_24pop_names

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
BN_all_rich <- colSums(allelic.richness(butternutgen_reorg)$Ar)/length(butternutgen_reorg@loc.n.all)	

##create data frame 
butternut_stat_df <- signif(cbind(butternut_mean_lon_lat[,2:3], BN_ind, BN_nall, BN_all_rich, BN_hexp),3)

##name columns and rows 
rownames(butternut_stat_df) <- butternut_24pop_names
colnames(butternut_stat_df) <- c("Mean Longitude", "Mean Latitude", "Number of Individuals", "MLG","Number of Alleles", "Allelic Richness", "Expected Heterozygosity")

##write out csv 
write.csv(butternut_stat_df, "genetic_analyses_results\\butternut_stat_df.csv")

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

