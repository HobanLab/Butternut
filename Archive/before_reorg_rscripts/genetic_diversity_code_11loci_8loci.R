
################################
########### Libraries ##########
################################

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
library(strataG)
library(plotrix)

###############################################
############## Working Directory ##############
###############################################

setwd("G:/My Drive/Hoban_Lab_Docs/Projects/Butternut_JUCI")

############################################################
######### Load Loci Documents - Genepop 11 loci ############
############################################################
##load in documents
butternutgen_11loci <- read.genepop("DataFiles/24Populations/11loci/butternut_24pops_11loci.gen", ncode = 3)

##reduce based on missing data 
butternutgen_nomd <- missingno(butternutgen_11loci, type = "geno", cutoff = 0.25, quiet = FALSE, freq = FALSE)

##data frame to reduce loading
butternut_latlon <- read.csv("DataFiles/24Populations/11loci/butternut_24pops_lonlat.csv")

##first reduce it by missing data 
butternut_latlon_nomd <- butternut_latlon[butternut_latlon$Ind %in% rownames(butternutgen_nomd@tab),]

##pop names
butternut_24pop_names <- unique(butternut_latlon_nomd$Pop)

###############################################
########## Loci Reduction - 8 Loci ############
###############################################
##remove loci - B157, B212_1, B262
butternutgen_8loci <- butternutgen_nomd[loc=-c(5,6,10)]


###############################################
############# Set up matrices #################
###############################################
##create data matrix
poppr_11loci <- poppr(butternutgen_nomd)
poppr_8loci <- poppr(butternutgen_8loci)

##allelic richness
allrich_11loci <- colSums(allelic.richness(butternutgen_nomd)$Ar)/length(butternutgen_nomd@loc.n.all)
allrich_8loci <- colSums(allelic.richness(butternutgen_8loci)$Ar)/length(butternutgen_nomd@loc.n.all)

##heterozygosity
hexp_11loci <- poppr_11loci[1:24,10]
hexp_8loci <- poppr_8loci[1:24,10]

##create data matrices
mean_coord_df <- read.csv("DataFiles\\24Populations\\mean_coords.csv")

##rename row names
rownames(mean_coord_df) <- mean_coord_df$X

##Reorganize document for moving 
mean_coord_df <- mean_coord_df[butternut_24pop_names,]
mean_coord_df <- mean_coord_df[,-1]

##heterozygosity
hexp_lat_11loci_df <- cbind(mean_coord_df$Mean.Lat, hexp_11loci)
hexp_lat_8loci_df <- cbind(mean_coord_df$Mean.Lat, hexp_8loci)

##allelic richness 
allrich_11loci_df <- cbind(mean_coord_df$Mean.Lat, allrich_11loci)
allrich_8loci_df <- cbind(mean_coord_df$Mean.Lat, allrich_8loci)

##name dfs
rownames(hexp_lat_11loci_df) <- butternut_24pop_names
rownames(hexp_lat_8loci_df) <- butternut_24pop_names
rownames(allrich_11loci_df) <- butternut_24pop_names
rownames(allrich_8loci_df) <- butternut_24pop_names

##reorder
hexp_lat_11loci_df <- hexp_lat_11loci_df[c("31","568","1014","7917",
                         "9101113a","9101113b","151","170","125147","126147",
                         "171188","APA12PACb","BV", "BVT","BVTGMVT2","CNF", "FKN","GMVT1", "MC",    
                         "OZ12","SFNF12", "WHWI", "WWI","VSH123"),]

hexp_lat_8loci_df <- hexp_lat_8loci_df[c("31","568","1014","7917",
                                           "9101113a","9101113b","151","170","125147","126147",
                                           "171188","APA12PACb","BV", "BVT","BVTGMVT2","CNF", "FKN","GMVT1", "MC",    
                                           "OZ12","SFNF12", "WHWI", "WWI","VSH123"),]

allrich_11loci_df <- allrich_11loci_df[c("31","568","1014","7917",
                                         "9101113a","9101113b","151","170","125147","126147",
                                         "171188","APA12PACb","BV", "BVT","BVTGMVT2","CNF", "FKN","GMVT1", "MC",    
                                         "OZ12","SFNF12", "WHWI", "WWI","VSH123"),]

allrich_8loci_df <- allrich_8loci_df[c("31","568","1014","7917",
                                         "9101113a","9101113b","151","170","125147","126147",
                                         "171188","APA12PACb","BV", "BVT","BVTGMVT2","CNF", "FKN","GMVT1", "MC",    
                                         "OZ12","SFNF12", "WHWI", "WWI","VSH123"),]

##create data frame
hexp_lat_11loci_df <- data.frame(hexp_lat_11loci_df)

##add color column
hexp_lat_11loci_df$Col <- NA

hexp_lat_11loci_df[1:6,3] <- "firebrick1"
hexp_lat_11loci_df[7:11,3] <- "firebrick4"
hexp_lat_11loci_df[12:24,3] <- "dodgerblue"

##
hexp_lat_8loci_df <- data.frame(hexp_lat_8loci_df)
hexp_lat_8loci_df$Col <- hexp_lat_11loci_df[,3]

allrich_11loci_df <- data.frame(allrich_11loci_df)
allrich_11loci_df$Col <- hexp_lat_11loci_df[,3]

allrich_8loci_df <- data.frame(allrich_8loci_df)
allrich_8loci_df$Col <- hexp_lat_11loci_df[,3]

###############################################
############# Write out Plots #################
###############################################
##now plot the comparison
hexp_lat_11loci_lm <- lm(hexp_lat_11loci_df[,2] ~ hexp_lat_11loci_df[,1])
hexp_lat_11loci_lm_sum <- summary(hexp_lat_11loci_lm)

##create r2 value 
hexp_11loci_r2 <- hexp_lat_11loci_lm_sum$adj.r.squared
hexp_11loci_coef <- hexp_lat_11loci_lm_sum$coefficients
hexp_11loci_pvalue <- hexp_11loci_coef[2,4]
hexp_rp_11loci <- vector('expression',2)
hexp_rp_11loci[1] <- substitute(expression(italic(R)^2 == MYVALUE), 
                              list(MYVALUE = format(hexp_11loci_r2,dig=3)))[2]
hexp_rp_11loci[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                             list(MYOTHERVALUE = format(hexp_11loci_pvalue, digits = 2)))[2]

##document write
pdf("Graphical_Stat_Results\\PostIndRemoval\\24pop\\11loci\\hexp_lat_11loci.pdf", width=8, height=6)

##write out plot
plot(hexp_lat_11loci_df[,1],hexp_lat_11loci_df[,2], xlab = "Mean Lat", ylab = "Expected Heterozygosity", pch = 16, cex = c(poppr_11loci[1:24,2]/75), col = hexp_lat_11loci_df[,3], ylim = c(0.75, 0.9))
legend('topleft', legend = hexp_rp_11loci, bty = 'n', border = "black", pt.cex = 1, cex = 0.8)
text(hexp_lat_11loci_df[,1],hexp_lat_11loci_df[,2], label = rownames(hexp_lat_11loci_df), pos = 2.5, offset = 0.5, cex = 0.8)
abline(hexp_lat_11loci_lm, col = "blue")
legend('top', legend = c("NB", "Ontario","US"), pch = 16, col = c("firebrick1","firebrick4","dodgerblue"))

dev.off()


##hexp 8 loci 
hexp_lat_8loci_lm <- lm(hexp_lat_8loci_df[,2] ~ hexp_lat_8loci_df[,1])
hexp_lat_8loci_lm_sum <- summary(hexp_lat_8loci_lm)

##create r2 value 
hexp_8loci_r2 <- hexp_lat_8loci_lm_sum$adj.r.squared
hexp_8loci_coef <- hexp_lat_8loci_lm_sum$coefficients
hexp_8loci_pvalue <- hexp_8loci_coef[2,4]
hexp_rp_8loci <- vector('expression',2)
hexp_rp_8loci[1] <- substitute(expression(italic(R)^2 == MYVALUE), 
                                list(MYVALUE = format(hexp_8loci_r2,dig=3)))[2]
hexp_rp_8loci[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                               list(MYOTHERVALUE = format(hexp_8loci_pvalue, digits = 2)))[2]

##document write
pdf("Graphical_Stat_Results\\PostIndRemoval\\24pop\\8loci\\hexp_lat_8loci.pdf", width=8, height=6)

##write out plot
plot(hexp_lat_8loci_df[,1],hexp_lat_8loci_df[,2], xlab = "Mean Lat", ylab = "Expected Heterozygosity 8 loci", pch = 16, cex = c(poppr_11loci[1:24,2]/75), col = hexp_lat_8loci_df[,3], ylim = c(0.75, 0.9))
legend('topleft', legend = hexp_rp_8loci, bty = 'n', border = "black", pt.cex = 1, cex = 0.8)
text(hexp_lat_8loci_df[,1],hexp_lat_8loci_df[,2], label = rownames(hexp_lat_8loci_df), pos = 2.5, offset = 0.5, cex = 0.8)
abline(hexp_lat_8loci_lm, col = "blue")
legend('top', legend = c("NB", "Ontario","US"), pch = 16, col = c("firebrick1","firebrick4","dodgerblue"))

dev.off()

##all rich 11 loci
##now plot the comparison
allrich_lat_11loci_lm <- lm(allrich_11loci_df[,2] ~ allrich_11loci_df[,1])
allrich_lat_11loci_lm_sum <- summary(allrich_lat_11loci_lm)

##create r2 value 
allrich_11loci_r2 <- allrich_lat_11loci_lm_sum$adj.r.squared
allrich_11loci_coef <- allrich_lat_11loci_lm_sum$coefficients
allrich_11loci_pvalue <- allrich_11loci_coef[2,4]
allrich_rp_11loci <- vector('expression',2)
allrich_rp_11loci[1] <- substitute(expression(italic(R)^2 == MYVALUE), 
                                list(MYVALUE = format(allrich_11loci_r2,dig=3)))[2]
allrich_rp_11loci[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                               list(MYOTHERVALUE = format(allrich_11loci_pvalue, digits = 2)))[2]

##document write
pdf("Graphical_Stat_Results\\PostIndRemoval\\24pop\\11loci\\allrich_lat_11loci.pdf", width=8, height=6)

##write out plot
plot(allrich_11loci_df[,1],allrich_11loci_df[,2], xlab = "Mean Lat", ylab = "Allelic Richness", main = "Allelic Richness 11 Loci",
     pch = 16, cex = c(poppr_11loci[1:24,2]/75), col = allrich_11loci_df[,3], ylim = c(6,10))
legend('topleft', legend = allrich_rp_11loci, bty = 'n', border = "black", pt.cex = 1, cex = 0.8)
text(allrich_11loci_df[,1],allrich_11loci_df[,2], label = rownames(allrich_11loci_df), pos = 2.5, offset = 0.5, cex = 0.8)
abline(allrich_lat_11loci_lm, col = "blue")
legend('bottomleft', legend = c("NB", "Ontario","US"), pch = 16, col = c("firebrick1","firebrick4","dodgerblue"))

dev.off()

##all rich 8 loci
##now plot the comparison
allrich_lat_8loci_lm <- lm(allrich_8loci_df[,2] ~ allrich_8loci_df[,1])
allrich_lat_8loci_lm_sum <- summary(allrich_lat_8loci_lm)

##create r2 value 
allrich_8loci_r2 <- allrich_lat_8loci_lm_sum$adj.r.squared
allrich_8loci_coef <- allrich_lat_8loci_lm_sum$coefficients
allrich_8loci_pvalue <- allrich_8loci_coef[2,4]
allrich_rp_8loci <- vector('expression',2)
allrich_rp_8loci[1] <- substitute(expression(italic(R)^2 == MYVALUE), 
                                   list(MYVALUE = format(allrich_8loci_r2,dig=3)))[2]
allrich_rp_8loci[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                                  list(MYOTHERVALUE = format(allrich_8loci_pvalue, digits = 2)))[2]

##document write
pdf("Graphical_Stat_Results\\PostIndRemoval\\24pop\\8loci\\allrich_lat_8loci.pdf", width=8, height=6)

##write out plot
plot(allrich_8loci_df[,1],allrich_8loci_df[,2], xlab = "Mean Lat", ylab = "Allelic Richness", main = "Allelic Richness 8 Loci",
     pch = 16, cex = c(poppr_11loci[1:24,2]/75), col = allrich_8loci_df[,3], ylim = c(4,7))
legend('topleft', legend = allrich_rp_8loci, bty = 'n', border = "black", pt.cex = 1, cex = 0.8)
text(allrich_8loci_df[,1],allrich_8loci_df[,2], label = rownames(allrich_8loci_df), pos = 2.5, offset = 0.5, cex = 0.8)
abline(allrich_lat_8loci_lm, col = "blue")
legend('bottomleft', legend = c("NB", "Ontario","US"), pch = 16, col = c("firebrick1","firebrick4","dodgerblue"))

dev.off()
N