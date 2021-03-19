######################################
########### Load Libraries ###########
######################################

library(adegenet)
library(poppr)
library(pegas)

#####################################
############# Load Files ############
#####################################
##set working directory
setwd("G:/My Drive/Hoban_Lab_Docs/Projects/Butternut_JUCI")

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

BN_poppr <- poppr(butternutgen_reorg)

##
BN_hexp <- BN_poppr[1:24, 10]
bn_hexp_df <- data.frame(butternut_24pop_names,BN_hexp)


##null alleles 
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
write.csv(bn_hwe, "G:\\Shared drives\\Emily_Schumacher\\Butternut\\butternut_publication_figures\\bn_hwe.csv")
write.csv(bn_hwe_reorg, "G:\\Shared drives\\Emily_Schumacher\\Butternut\\butternut_publication_figures\\bn_hwe_reorg.csv")
write.csv(ld_comp_df, "G:\\Shared drives\\Emily_Schumacher\\Butternut\\butternut_publication_figures\\ld_loci.csv")

