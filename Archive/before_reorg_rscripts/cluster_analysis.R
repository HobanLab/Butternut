######################
###### Libraries #####
######################

library(adegenet)

###############################################
################# Load in Files  ##############
###############################################
##set working directory
setwd("C:\\Users\\eschumacher\\Documents\\GitHub\\butternut")

###############################################
################ run PCoAs ####################
###############################################
##PCoA list 
related_red_list <- list.files(path = "DataFiles\\before_reorg", pattern = "_rel_red.csv$")

##list out 
relatedness_reduced_df_list <- list()

##genind list 
gen_list <- list.files(path = "DataFiles\\before_reorg", pattern = "24pops.gen$")

##create genind list 
genind_list <- list()

##reduce genind
genind_red_list <- list()

##convert to genpop
genpop_list <- list()

##create list for the tab documents
pco_tab_list <- list()

##pcoa list
pco_list <- list()

##loop to create PCoA plots
for(i in 1:length(related_red_list)){
  
  ##read in csv files 
  relatedness_reduced_df_list[[i]] <- read.csv(paste0("DataFiles\\before_reorg\\", related_red_list[[i]]))
  
}

##24 population names 
butternut_24pop_names <- unique(relatedness_reduced_df_list[[1]]$Pop)

##load in lon/lat data frame 
butternut_lonlat_8loci_nomd <- read.csv("DataFiles\\before_reorg\\butternut_lonlat_8loci_nomd.csv")

##reduce genind files 
for(i in 1:length(gen_list)){
  
  ##read in genind files
  genind_list[[i]] <- read.genepop(paste0("DataFiles\\before_reorg\\", gen_list[[i]]), ncode = 3)
  
  
}

##name rows in the 8 loci document 
rownames(genind_list[[2]]@tab) <- butternut_lonlat_8loci_nomd$Pop

##reduce genind files 
for(i in 1:length(gen_list)){
  
  ##read in genind files
  genind_red_list[[i]] <- genind_list[[i]][!rownames(genind_list[[i]]@tab) %in% relatedness_reduced_df_list[[i]][,1] ,]
  
  ##create genpop files
  genpop_list[[i]] <- genind2genpop(genind_red_list[[i]])
  
  ##prepare doc 
  pco_tab_list[[i]] <- tab(genpop_list[[i]], freq = TRUE, NA.method = "mean")
  
  ##run PCoA
  pco_list[[i]] <- dudi.pco(dist(pco_tab_list[[i]]), scannf = FALSE, nf = 2)
  
}

##Calculate percent variation explained 
sum_eig <- sum(butternut_reorg_pco$eig)
pc1 <- (butternut_reorg_pco$eig[[1]]/sum_eig)*100
pc2 <- (butternut_reorg_pco$eig[[2]]/sum_eig)*100

##pcoa data frame
pco_df <- data.frame(cbind(pco1$li$A1, pco1$li$A2))
pco_df$col <- NA
pco_df[1:504,3] <- "firebrick1"
pco_df[505:977,3] <- "dodgerblue"

butternut_pco_nb <- rbind(butternut_reorg_pco$li[c("31","568","1014","7917","9101113a","9101113b"),1:2])
butternut_pco_ot <- rbind(butternut_reorg_pco$li[c("151","170","125147","126147","171188"),1:2])

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