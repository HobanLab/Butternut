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