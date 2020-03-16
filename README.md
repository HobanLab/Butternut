# butternut
Collection of all code for butternut
library(adegenet)
##load butternut document

##setwd("C:/Users/eschumacher/Documents/MortonArboretum/")
##personal computer
setwd("C:/Users/eksch/Documents/MortonArboretum")

butternutgen <- read.genepop("butternut_USCAN.gen", ncode = 3)

butternutpop <- genind2genpop(butternutgen)

##first identify if there were scoring discrepencies between Sean and Jeanne

loci <- c("B114","B159","WGA","A5_2","B157","B212_2","B121",	"B147",	"B249",	"B262","B264")

pdf(file=paste0("locibarplot.pdf"),width=40,height=9)

for(a in loci){
  
  new <- butternutpop[,which(grepl(a,colnames(butternutpop@tab)))]@tab
  
  barplot(new, las = 2, beside = TRUE, col = c("red", "blue"), legend.text =  c("Canada", "United States"))
  
}
dev.off()
