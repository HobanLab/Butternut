#This code removes individuals based on missing data and relatedness. 
#This data file starts with the 1,721 individuals following geographic 
#reorganization into 24 populations and rebinning analysis.
##This generates the data file used for all diversity analyses
#Then, individuals were removed based on relatedness. 
#Relatedness between individuals can over-represent genotypes in clustering
#analysis and bias these analyses away from finding all the groupings. 
#Therefore we were able to use this code to remove individuals based on 
#greater than 25% similar genotypes using the Loiselle statistic
#before running PCoA and structure analysis.

##########################
######## Libraries #######
##########################

library(adegenet)
library(Demerelate)
library(poppr)

#####################################
############ Directories ############
#####################################
butternut_drive <- "C:\\Users\\eschumacher\\Documents\\GitHub\\butternut"

#####################################
############# Load Files ############
#####################################
setwd(butternut_drive)

##load current working reorganized genepop file 
butternutgen_reorg <- read.genepop("data_files\\after_reorg\\postgeo_rebin_24pop.gen", ncode = 3)

##reduce genind file for individuals with greater than 25% missing data 
butternutgen_nomd <- missingno(butternutgen_reorg, type = "geno", cutoff = 0.25, quiet = FALSE, freq = FALSE)

##load relatedness file 
reorg_relatedness <- read.csv("data_files\\after_reorg\\reorg_relatedness.csv")

###################################################################
############# Reduced Individuals based on Relatedness ############
###################################################################
##rename individuals in genind file 
rownames(butternutgen_nomd@tab) <- reorg_relatedness$Ind

##create pop name code 
butternut_24pop_names <- unique(reorg_relatedness$Pop)

##rename populations in genind no missing data file 
levels(butternutgen_nomd@pop) <- butternut_24pop_names

####reduce relatedness
##run relatedness code
reorg_relate_df <- Demerelate(reorg_relatedness, object = T, value = "loiselle")

##now identify how many individuals have greater than 25% relatedness = half siblings
butternut_halfsib_names <- names(which(unlist(reorg_relate_df$Empirical_Relatedness) > 0.25))

##then use this to create a document which has all the unique individual numbers for every highly related individuals
butternut_halfsib_names_cleanfront <- gsub("^.*\\.","", butternut_halfsib_names)

butternut_halfsib_names_cleanback <- gsub("^.*\\_","", butternut_halfsib_names_cleanfront)

relate_ind_remove <- unique(butternut_halfsib_names_cleanback)

##then subset genind file
butternutgen_relatedness_reduced <- butternutgen_nomd[!rownames(butternutgen_nomd@tab) %in% relate_ind_remove,]

##subset data frame 
reorg_relatedness_reduced <- reorg_relatedness[!reorg_relatedness$Ind %in% relate_ind_remove,]

###name pops
levels(butternutgen_relatedness_reduced@pop) <- butternut_24pop_names
 
##write out genind file in genalex format
setwd("data_files\\after_reorg")
butternutgen_relatedness_reduced <- genind2genalex(butternutgen_relatedness_reduced, 
                                                   filename = "butternutgen_relatedness_reduced.csv")
write.csv(reorg_relatedness_reduced, "reorg_relatedness_reduced.csv")

##write out no missing data genind file and then convert to genind in genalex
butternutgen_nomd <- genind2genalex(butternutgen_nomd, filename = "butternutgen_nomd.csv")
