##########################
######## Libraries #######
##########################

library(diveRsity)
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
library(plotrix)
library(ggpmisc)
library(factoextra)
library(GISTools)
library(raster)
library(rgdal)
library(sp)
library(PopGenReport)

#####################################
############ Directories ############
#####################################
shared_drive <- "G:\\Shared drives\\Emily_Schumacher\\butternut_publication_figures"
butternut_drive <- "G:\\My Drive\\Hoban_Lab_Docs\\Projects\\Butternut_JUCI"

#####################################
############ Null Alleles ###########
#####################################
##setwd
setwd(butternut_drive)

##load current working reorganized genepop file 
butternutgen_reorg <- read.genepop("DataFiles\\24Populations\\reorg\\reorg_gen_24pop.gen", ncode = 3)

##Calculate null alleles
null_all <- null.all(butternutgen_reorg)

##butternut loci
loci_names <- c("B114","B159","WGA","A5_2", "B157","B212_2",
                "B121","B147","B249", "B262", "B264")  
colnames(null_all$null.allele.freq$summary1) <- loci_names

write.csv(null_all$null.allele.freq$summary1, "Null_alleles_Freq.csv")

###Null allele calculations 
pdf("Null_Alleles_by_Loci.pdf", width = 10, height = 8)
barplot(null_all$null.allele.freq$summary1[2,], xlab = "Loci", ylab = "Median Percentage of Null Alleles", main = "Null Alleles by Loci", ylim = c(0,1), col = "dodgerblue")
dev.off()
