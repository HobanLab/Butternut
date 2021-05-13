#This code is used to generate PCoAs for the analysis
#This file uses the reduced individual document from the 
#relatedness removal document 

##########################
######## Libraries #######
##########################

library(adegenet)
library(diveRsity)

#####################################
############ Directories ############
#####################################
##set directory to all butternut files 
butternut_drive <- "C:\\Users\\eschumacher\\Documents\\GitHub\\butternut"

#####################################
############# Load Files ############
#####################################
setwd(butternut_drive)

##load current reorganized and reduced genepop with reduced by relatedness  
#only use if not converted 
#butternutgen_red <- arp2gen("data_files\\after_reorg\\butternutgen_relatedness_reduced.arp")

##load genind file 
butternutgen_relate <- read.genepop("data_files\\after_reorg\\butternutgen_relatedness_reduced.gen", ncode = 3)

##load relatedness data frame 
reorg_relatedness_df <- read.csv("data_files\\after_reorg\\reorg_relatedness_reduced.csv")

##create pop name code 
butternut_24pop_names <- unique(reorg_relatedness_df$Pop)

##############################
######### PCOA Code ##########
##############################
##now reorg into a genepop file 
butternutpop_relate <- genind2genpop(butternutgen_relate)

##run the PCA function
butternut_reorg_pco <- dudi.pco(dist.genpop(butternutpop_relate, meth = 2), nf = 2, scannf = FALSE)

##format df
rownames(butternut_reorg_pco$li) <- butternut_24pop_names

##re-org data for color coding
butternut_pco_nb <- rbind(butternut_reorg_pco$li[c("31","568","1014","7917","9101113a","9101113b"),1:2])
butternut_pco_ot <- rbind(butternut_reorg_pco$li[c("151","125147","126147"),1:2])
butternut_pco_qu <- rbind(butternut_reorg_pco$li[c("170","171188"),1:2])

##Calculate percent variation explained 
sum_eig <- sum(butternut_reorg_pco$eig)
pc1 <- (butternut_reorg_pco$eig[[1]]/sum_eig)*100
pc2 <- (butternut_reorg_pco$eig[[2]]/sum_eig)*100

##PCOA of the reorg
pdf("genetic_analyses_results\\PCoA.pdf", width = 8, height = 6)
##plot New Brunswick populations
plot(butternut_pco_nb$A1, butternut_pco_nb$A2, pch = 17, 
     xlab = paste0("PC1", sep = " ", "(",round(pc1, digits = 1), "%", ")"), 
     ylab = paste0("PC2", sep = " ", "(",round(pc2, digits = 1), "%",")"), 
     main = "PCoA 24 Populations Reorg", col = "firebrick1", xlim = c(-0.3, 0.3), ylim = c(-0.3, 0.3))
##name with populations 
text(butternut_pco_nb$A1, butternut_pco_nb$A2, label = rownames(butternut_pco_nb), pos = 2, cex = 0.8)
##plot Ontario populations
points(butternut_pco_ot$A1, butternut_pco_ot$A2,pch = 17, col = "firebrick4")
text(butternut_pco_ot$A1, butternut_pco_ot$A2, label = rownames(butternut_pco_ot), pos = 3, cex = 0.8)
##plot Quebec populations 
points(butternut_pco_qu$A1, butternut_pco_qu$A2,pch = 17, col = "lightsalmon")
text(butternut_pco_qu$A1, butternut_pco_qu$A2, label = rownames(butternut_pco_qu), pos = 3, cex = 0.8)
##plot US populations 
points(butternut_reorg_pco[12:24,]$li$A1, butternut_reorg_pco[12:24,]$li$A2, col = "dodgerblue", pch = 17)
text(butternut_reorg_pco[12:24,]$li$A1, butternut_reorg_pco[12:24,]$li$A2, 
     label = rownames(butternut_reorg_pco[12:24,]$li), pos = 4, cex = 0.8)
##create legend 
legend('topleft', pch = 17, col = c("firebrick1","firebrick4", "lightsalmon", "dodgerblue"), 
       legend = c("New Brunswick", "Ontario", "Quebec", "United States"))
abline(h = 0)
abline(v = 0)
dev.off()
