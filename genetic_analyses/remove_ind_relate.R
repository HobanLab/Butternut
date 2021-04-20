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
butternutgen_reorg <- read.genepop("data_files\\after_reorg\\reorg_gen_24pop.gen", ncode = 3)

##load relatedness file 
reorg_relatedness <- read.csv("data_files\\after_reorg\\reorg_relatedness.csv")

###################################################################
############# Reduced Individuals based on Relatedness ############
###################################################################

##rename individuals in genind file 
rownames(butternutgen_reorg@tab) <- reorg_relatedness$Ind

##create pop name code 
butternut_24pop_names <- unique(reorg_relatedness$Pop)

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
butternutgen_relatedness_reduced <- butternutgen_reorg[!rownames(butternutgen_reorg@tab) %in% relate_ind_remove,]

###name pops
levels(butternutgen_relatedness_reduced@pop) <- butternut_24pop_names
 
##write out genind file in genalex format
setwd("data_files\\after_reorg")
butternutgen_relatedness_reduced <- genind2genalex(butternutgen_relatedness_reduced, 
                                                   filename = "butternutgen_relatedness_reduced.csv")

