##This script is used to determine if there are differences between 

#####################
#     Libraries     #
#####################

library(adegenet)

###############################
#     Scoring Comparison      #
###############################
#set working directory
setwd("../data_files/before_reorg")

#three pop genind and genpop files
butternut_3pop_gen <- read.genepop("butternut_3pops.gen", ncode = 3)
butternut_3pop_pop <- genind2genpop(butternut_3pop_gen)

#create list of loci 
loci <- c("B114","B159","WGA","A5_2","B157","B212_2","B121",	"B147",	"B249",	"B262","B264")

#loop to compare scoring between Jeanne and Sean
pdf("butternut_3pop_scoring_barplot.pdf",width=40,height=9)

for(a in loci){
  
  #separate by loci
  butternut_3pop_scoring <- butternut_3pop_pop[,which(grepl(a,colnames(butternut_3pop_pop@tab)))]@tab
  
  #calculate scoring frequency for each allele by researcher scoring
  for(p in 1:3) butternut_3pop_scoring[p,]<-butternut_3pop_scoring[p,]/sum(butternut_3pop_scoring[p,])
  
  #plot frequency for each allele
  butternut_3pop_scoring_barplot <- barplot(butternut_3pop_scoring, las = 2, beside = TRUE, col = c("firebrick4","firebrick2","dodgerblue"), legend.text =  c("Canada_Sean", "Canada_Jean","United States"), ylim = c(0,1), main = paste0(a))
  
}

dev.off()

sessionInfo()