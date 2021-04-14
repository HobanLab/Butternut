######################
###### Libraries #####
######################

library(adegenet)
library(tidyr)
library(poppr)

###############################################
############## Working Directory ##############
###############################################

setwd("C:\\Users\\eschumacher\\Documents\\GitHub\\butternut")

###############################################
############# Individual Reduction ############
###############################################
##following pop reclassification - load 24 pop document
butternutgen_11loci <- read.genepop("DataFiles\\butternut_24pops_11loci.gen", ncode = 3)

##now remove individuals with greater than 25% missing data
butternutgen_nomd_11loci <- missingno(butternutgen_11loci, type = "geno", cutoff = 0.25, quiet = FALSE, freq = FALSE)

#load in relatedness file
butternut_rel_11loci <- read.csv("")

##calculate relatedness
butternutrel <- Demerelate(butternut_24pop_rel, object = T, value = "loiselle")

##now identify how many individuals have greater than 25% relatedness = half siblings
butternut_halfsib_names <- names(which(unlist(butternutrel$Empirical_Relatedness) > 0.25))

##then use this to create a document which has all the unique individual numbers for every highly related individuals
butternut_halfsib_names_cleanfront <- gsub("^.*\\.","", butternut_halfsib_names)

butternut_halfsib_names_cleanback <- gsub("^.*\\_","", butternut_halfsib_names_cleanfront)

relate_ind_remove <- unique(butternut_halfsib_names_cleanback)

##then subset genind file
butternutgen_24pop_reduced <- butternutgen_nomd[!rownames(butternutgen_nomd@tab) %in% relate_ind_remove,]

##data frame to reduce loading
butternut_latlon <- read.csv("G:/My Drive/Hoban_Lab_Docs/Projects/Butternut_JUCI/DataFiles/24Populations/11loci/butternut_24pops_lonlat.csv")

##first reduce it by missing data 
butternut_latlon_nomd <- butternut_latlon[butternut_latlon$Ind %in% butternut_24pop_rel$Ind,]

##then reduce it by relatedness individuals
butternut_24pop_coords <- butternut_latlon_nomd[!butternut_latlon_nomd$Ind %in% relate_ind_remove,]

##create a genpop file
butternutpop <- genind2genpop(butternutgen_24pop_reduced)

###############################################
################ run PCoAs ####################
###############################################

##24 pops
butternut_24pop_names <- unique(butternut_latlon$Pop)

##prepare doc 
pco_prep <- tab(butternutgen_24pop_reduced, freq = TRUE, NA.method = "mean")

##run PCoA
pco1 <- dudi.pco(dist(pco_prep), scannf = FALSE, nf = 2)

##pcoa data frame
pco_df <- data.frame(cbind(pco1$li$A1, pco1$li$A2))
pco_df$col <- NA
pco_df[1:504,3] <- "firebrick1"
pco_df[505:977,3] <- "dodgerblue"

##now plot 
pdf("Graphical_Stat_Results\\PostIndRemoval\\24pop\\indlevel_PCoA_10_14_20.pdf")
plot(pco_df[,1], pco_df[,2], col = pco_df[,3], pch = 17, ylab = "PC2 (16.2%)", xlab = "PC1 (20.2%)")
legend("topright", col = c("firebrick1","dodgerblue"), pch = 17, legend = c("Canada", "US"))
abline(h = 0, col = "black")
abline(v = 0, col = "black")
dev.off()

##now do population level
pco_pop_prep <- tab(butternutpop, freq = TRUE, NA.method = "mean")

##run pcoa on genpop file
pco_pop <- dudi.pco(dist(pco_pop_prep), scannf = FALSE, nf = 2)

##create a data frame
pco_pop_df <- data.frame(cbind(pco_pop$li$A1, pco_pop$li$A2))
rownames(pco_pop_df) <- butternut_24pop_names

pco_pop_df <- pco_pop_df[c("31","568","1014","7917",
   "9101113a","9101113b","151","170","125147","126147",
   "171188","APA12PACb","BV", "BVT","BVTGMVT2","CNF", "FKN","GMVT1", "MC",    
   "OZ12","SFNF12", "WHWI", "WWI","VSH123"),]

pco_pop_df$Col <- NA

pco_pop_df[1:6,3] <- "firebrick1"
pco_pop_df[7:11,3] <- "firebrick4"
pco_pop_df[12:24,3] <- "dodgerblue"

##write out PDF
pdf("Graphical_Stat_Results\\PostIndRemoval\\24pop\\POPlevel_PCoA_10_14_20.pdf")
plot(pco_pop_df[,1], pco_pop_df[,2], col = pco_pop_df[,3], pch = 17, xlab = "PC1 (11.6%)", ylab = "PC2 (6.03%)",
     main = "Pop Level PCoA")
text(pco_pop_df[,1], pco_pop_df[,2], label = rownames(pco_pop_df), cex = 0.5, pos = 2)
legend("topright", col = c("firebrick1","firebrick4","dodgerblue"), pch = 17, legend = c("NB", "Ontario","US"))
abline(h = 0)
abline(v = 0)
dev.off()