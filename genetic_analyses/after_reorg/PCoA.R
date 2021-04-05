##########################
######## Libraries #######
##########################

library(adegenet)
library(stringr)
library(tidyr)
library(poppr)
library(Demerelate)
library(PopGenReport)

#####################################
############ Directories ############
#####################################
butternut_drive <- "G:\\My Drive\\Hoban_Lab_Docs\\Projects\\Butternut_JUCI"

#####################################
############# Load Files ############
#####################################
setwd(butternut_drive)

################################################################
############### Genetic Structure Analyses #####################
################################################################
##load current working reorganized genepop file 
butternutgen_reorg <- read.genepop("DataFiles\\24Populations\\reorg\\reorg_gen_24pop.gen", ncode = 3)

##load relatedness file 
reorg_relatedness <- read.csv("DataFiles\\24Populations\\reorg\\reorg_relatedness.csv")

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

###################Begin PCoA code
butternutpop_reorg <- genind2genpop(butternutgen_relatedness_reduced)

##run the PCA function
butternut_reorg_pco <- dudi.pco(dist.genpop(butternutpop_reorg, meth = 2), nf = 2, scannf = FALSE)

##format df
rownames(butternut_reorg_pco$li) <- butternut_24pop_names

##re-org data for color coding
butternut_pco_nb <- rbind(butternut_reorg_pco$li[c("31","568","1014","7917","9101113a","9101113b"),1:2])
butternut_pco_ot <- rbind(butternut_reorg_pco$li[c("151","170","125147","126147","171188"),1:2])

##Calculate percent variation explained 
sum_eig <- sum(butternut_reorg_pco$eig)
pc1 <- (butternut_reorg_pco$eig[[1]]/sum_eig)*100
pc2 <- (butternut_reorg_pco$eig[[2]]/sum_eig)*100

##PCOA of the reorg
pdf("Graphical_Stat_Results\\PostIndRemoval\\24pop\\Reorg_Results\\STR\\PCoA.pdf", width = 8, height = 8)
##plot New Brunswick populations
plot(butternut_pco_nb$A1, butternut_pco_nb$A2, pch = 17, 
     xlab = paste0("PC1", sep = " ", "(",round(pc1, digits = 1), "%", ")"), 
     ylab = paste0("PC2", sep = " ", "(",round(pc2, digits = 1), "%",")"), 
     main = "PCoA 24 Populations Reorg", col = "firebrick1", xlim = c(-0.3, 0.3), ylim = c(-0.3, 0.3))
text(butternut_pco_nb$A1, butternut_pco_nb$A2, label = rownames(butternut_pco_nb), pos = 2, cex = 0.8)
##plot Ontario populations
points(butternut_pco_ot$A1, butternut_pco_ot$A2,pch = 17, col = "firebrick4")
text(butternut_pco_ot$A1, butternut_pco_ot$A2, label = rownames(butternut_pco_ot), pos = 3, cex = 0.8)
##plot US populations 
points(butternut_reorg_pco[12:24,]$li$A1, butternut_reorg_pco[12:24,]$li$A2, col = "dodgerblue", pch = 17)
text(butternut_reorg_pco[12:24,]$li$A1, butternut_reorg_pco[12:24,]$li$A2, 
     label = rownames(butternut_reorg_pco[12:24,]$li), pos = 4, cex = 0.8)
##create legend 
legend('topleft', pch = 17, col = c("firebrick1","firebrick4", "dodgerblue"), 
       legend = c("New Brunswick", "Ontario", "United States"))
abline(h = 0)
abline(v = 0)
dev.off()
