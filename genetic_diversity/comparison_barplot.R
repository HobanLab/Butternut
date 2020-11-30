##########################
######## Libraries #######
##########################

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

###############################################
############# Scoring Comparison ##############
###############################################

##three pop genind and genpop files
butternut_3pop_gen <- read.genepop("G:/My Drive/Hoban_Lab_Docs/Projects/Butternut_JUCI/DataFiles/butternut3.gen", ncode = 3)
butternut_3pop_pop <- genind2genpop(butternut_3pop_gen)

##create list of loci 
loci <- c("B114","B159","WGA","A5_2","B157","B212_2","B121",	"B147",	"B249",	"B262","B264")

##loop to compare scoring between Jeanne and Sean
pdf(file=paste0("butternut_3pop_scoring_barplot.pdf"),width=40,height=9)

for(a in loci){
  
  butternut_3pop_scoring <- butternut_3pop_pop[,which(grepl(a,colnames(butternut_3pop_pop@tab)))]@tab
  
  for(p in 1:3) butternut_3pop_scoring[p,]<-butternut_3pop_scoring[p,]/sum(butternut_3pop_scoring[p,])
  
  butternut_3pop_scoring_barplot <- barplot(butternut_3pop_scoring, las = 2, beside = TRUE, col = c("firebrick4","firebrick2","dodgerblue"), legend.text =  c("Canada_Sean", "Canada_Jean","United States"), ylim = c(0,1), main = paste0(a))
  
}

dev.off()