######################
###### Libraries #####
######################

library(adegenet)
library(tidyr)
library(poppr)
library(Demerelate)

###############################################
############## Working Directory ##############
###############################################

setwd("C:\\Users\\eschumacher\\Documents\\GitHub\\butternut")

###############################################
############# Individual Reduction ############
###############################################
##list of genind files for individual reduction
gen_list <- list.files(path = "DataFiles\\before_reorg", pattern = "_24pops.gen$")

##relatedness data frame list 
related_list <- list.files(path = "DataFiles\\before_reorg", pattern = "_rel.csv$")

##relatedness data frames list 
related_df_list <- list()

##create lists 
genind_list <- list()

##remove missing data 
genind_nomd_list <- list()

##relatedness 
demerelate_list <- list()

##half sibling names list 
halfsib_names_list <- list()

##create list for removing everything before the period 
cleanfront_list <- list()

##create list for removing 2nd individual 
cleanback_list <- list()

##create list for related individuals 
relate_ind_list <- list()

##reduce by relatedness
genind_rel_red <- list()

##load in lon lat documents 
lonlat_list <- list.files(path = "DataFiles\\before_reorg", pattern = "_lonlat.csv$")

##lonlat list to read in csv files 
lonlat_df_list <- list()

##lon lat reduction from missing data 
lonlat_nomd_list <- list()

##

##loop to read in genind files, convert and reduce 
for(i in 1:length(gen_list)){
  
  ##read in genind files 
  genind_list[[i]] <- read.genepop(paste0("DataFiles\\before_reorg\\", gen_list[[i]]), ncode = 3)
  
  ##remove individuals with greater than 25% missing data
  genind_nomd_list[[i]] <- missingno(genind_list[[i]], type = "geno", cutoff = 0.25, quiet = FALSE, freq = FALSE)
  
  ##read in csv of relatedness files
  related_df_list[[i]] <- read.csv(paste0("DataFiles\\before_reorg\\", related_list[[i]]))

  ##calculate relatedness
  demerelate_list[[i]] <- Demerelate(related_df_list[[i]], object = T, value = "loiselle")
  
  ##half sibling names from the demerelate
  halfsib_names_list[[i]] <- names(which(unlist(demerelate_list[[i]][2]) > 0.25))
  
  ##remove everything before the period 
  cleanfront_list[[i]] <- gsub("^.*\\.","", halfsib_names_list[[i]])
  
  ##remove 2nd individuals in comparison
  cleanback_list[[i]] <- gsub("^.*\\_","", cleanfront_list[[i]])
  
  ##lists for individuals to reduce by name  
  relate_ind_list[[i]] <- unique(cleanback_list[[i]])
  
  ##create row names 
  rownames(genind_nomd_list[[i]]@tab) <- related_df_list[[2]][,1]
 
  ##subset genind files 
  genind_rel_red[[i]] <- genind_nomd_list[[i]][!rownames(genind_nomd_list[[i]]@tab) %in% relate_ind_list[[i]],]
  
  ##load in lon/lat documents 
  lonlat_df_list[[i]] <- read.csv(paste0(paste0("DataFiles\\before_reorg\\",lonlat_list[[i]])))
  
  
  
}

##data frame to reduce loading
butternut_lonlat_11loci <- read.csv("DataFiles\\butternut_24pops_lonlat_11loci.csv")

##reduce lon lat data frame by missing data 
butternut_lonlat_nomd_11loci <- butternut_lonlat_11loci[as.character(butternut_lonlat_11loci$Ind ) %in% 
                                                          rownames(butternutgen_nomd_11loci@tab),]

##then reduce it by relatedness individuals
butternut_lonlat_rel_red_11loci <- butternut_lonlat_nomd_11loci[!butternut_lonlat_nomd_11loci$Ind %in% 
                                                                  relate_ind_remove,]

#############################################
############# Write out Files ###############
#############################################
setwd("genetic_analyses\\genetic_analyses_results")
###write out data frames
write.csv(butternut_red_11loci, "DataFiles\\before_reorg\\butternut_11loci_rel.csv")
write.csv(butternut_lonlat_nomd_11loci, "butternut_lonlat_nomd_11loci.csv")
write.csv(butternut_lonlat_rel_red_11loci, "butternut_lonlat_rel_red_11loci.csv")


