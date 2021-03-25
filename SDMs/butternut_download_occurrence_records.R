### DESCRIPTION:
###############This is based on code provide by Emily Beckman, RA at the Morton Arboretum (Github repository: https://github.com/esbeckman/IMLS_Beckman)
# This script provides instructions and code chunks for downloading wild
#   occurrence points from:
# GLOBAL DATABASES (though all likely have U.S. bias?)
# Global Biodiversity Information Facility (GBIF)
# Integrated Digitized Biocollections (iDigBio)
# U.S. Herbarium Consortia (SERNEC, SEINet, etc.)
# Botanical Information and Ecology Network (BIEN)
# NATIONAL DATABASES
# Forest Inventory and Analysis (FIA) Program of the USDA Forest Service

### INPUT:
# target_taxa_with_syn.csv (list of target taxa)
# columns:
# 1. "taxon_name" (genus, species, infra rank, and infra name, all
#    separated by one space each; hybrid symbol should be " x ", rather
#    than "_" or "???", and go between genus and species)
# 2. (optional) other data you want to keep with taxa info

### OUTPUTS:
# gbif_raw.csv
# idigbio_raw.csv
# herbaria_raw.csv
# bien_raw.csv
# fia_raw.csv

#################
### LIBRARIES ###
#################

library(plyr)
library(tidyverse) #ggplot2,dplyr,tidyr,readr,purrr,tibble,stringr,forcats
library(spocc)
library(rgbif)
library(data.table)
library(BIEN)
library(ridigbio)
library(batchtools)
library(googledrive)

################################################################################
# A) Read in target taxa list
################################################################################

setwd("G:/My Drive/Hoban_Lab_Docs/Projects/Butternut_JUCI/DataFiles/Occurrence")


# list of target taxon names
taxon_names <- c("Juglans cinerea")
# list of target species names
species_names <- c("Juglans cinerea")


################################################################################
# B) Global Biodiversity Information Facility (GBIF) download
################################################################################

# GBIF account user information
# if you don't have account yet, go to https://www.gbif.org then click
#   "Login" in top right corner, then click "Register"
# !!! FILL THIS IN WITH YOUR INFO:
user <- "ekschu"
pwd <- 
email <- 

# get GBIF taxon keys for all taxa in target list
keys <- sapply(taxon_names,function(x) name_backbone(name=x)$speciesKey,
               simplify = "array")

# remove duplicate and NULL keys
keys_nodup <- keys[!duplicated(keys) & keys != "NULL"]
# create data frame of keys and matching taxon_name

gbif_codes <- map_df(keys_nodup,~as.data.frame(.x),.id="taxon_name")
names(gbif_codes)[2] <- "speciesKey"
# create vector of keys as input into gbif download
gbif_taxon_keys <- gbif_codes[,2]

# download GBIF data (Darwin Core Archive format)
gbif_download <- occ_download(
  pred_in("taxonKey", gbif_taxon_keys),
  pred_in("basisOfRecord", c("PRESERVED_SPECIMEN",
                             "HUMAN_OBSERVATION","FOSSIL_SPECIMEN","OBSERVATION",
                             "UNKNOWN","MACHINE_OBSERVATION","MATERIAL_SAMPLE",
                             "LITERATURE")),
  #pred("hasCoordinate", TRUE),
  #pred("hasGeospatialIssue", FALSE),
  format = "DWCA", #"SIMPLE_CSV"
  user=user,pwd=pwd,
  email=email)
# load gbif data just downloaded
# create new folder for data and set as working directory
dir.create(file.path(getwd(),"gbif_read_in"))
setwd(file.path(getwd(),"gbif_read_in"))

# must wait for download to complete before continuing;
# it may take a while (up to 3 hours) if you have a large taxa list;
# you can check download status here: https://www.gbif.org/user/download

# download and unzip before reading in
gbif_download # !!! PASTE "Download key" as first argument in next two lines !!!
occ_download_get(key="0039229-200221144449610", overwrite=TRUE)
unzip("0039229-200221144449610.zip")
# read in data
gbif_raw <- fread("occurrence.txt",quote="")
nrow(gbif_raw) #2399719
# write file
setwd("./..")
write.csv(gbif_raw, "gbif_raw.csv")

##gbif 

gbif_raw <- read.csv("gbif_raw2.csv")

################################################################################
# C) Integrated Digitized Biocollections (iDigBio) download
################################################################################

# I'm not sure the R interface actually gets all fields available; use manual
#    way down below (starts line 144) if you want to be sure

# download iDigBio data for target taxa
# we have to go taxon by taxon; function can only return 100,000
#   records at once and Quercus has more than that so can't download by genera
idigbio_raw <- data.frame()
for(i in 1:length(taxon_names)){
  output_new <- idig_search_records(rq=list(scientificname=taxon_names[[i]]),
                                    fields="all")
  idigbio_raw <- rbind(idigbio_raw,output_new)
  print(paste(round(i/length(taxon_names)*100,digits=1),"% complete",sep=""))
}
nrow(idigbio_raw) #153467
# remove rows that are lists
idigbio_raw <- idigbio_raw %>% select(everything(),-commonnames,-flags,
                                      -mediarecords,-recordids)
# write file
write.csv(idigbio_raw, "idigbio_raw.csv")

##
idigbio_raw <- read.csv("idigbio_raw.csv")


################################################################################
# D) U.S. Herbaria Consortia (SERNEC, SEINet, etc.) download
################################################################################

# First, download raw data
# Go to http://sernecportal.org/portal/collections/harvestparams.php
# Type your target genus name into the "scientific name" box and click
#   "List Display"
# Click the Download Specimen Data button (arrow pointing down into a box),
#   in the top right corner
# In the pop-up window, select the "Darwin Core" radio button,
#   uncheck everything in the "Data Extensions:" section, and
#   select the "UTF-8 (unicode)" radio button
#   leave other fields as-is
# Click "Download Data"

# If you have more than one target genus, repeat the above steps for the
#   other genera

# Move all the zipped files you downloaded into a "sernec_read_in" folder
#   within your working directory
# Unzip each file and pull the "occurrences.csv" file out into the
#   "sernec_read_in" folder -- obviously "keep both" if prompted

# read in raw occurrence points

sernec_raw <- read.csv("G:/My Drive/Hoban_Lab_Docs/Projects/Butternut_JUCI/DataFiles/Occurrence/sernec_read_in/occurrences.csv")

sernec_raw <- read.csv("G:/My Drive/Hoban_Lab_Docs/Projects/Butternut_JUCI/DataFiles/Occurrence/sernec_read_in/occurrences.csv")


################################################################################
# E) Botanical Information and Ecology Network (BIEN) download
################################################################################

# information about functions in package
#vignette("BIEN")

# download BIEN occurrence data for target genera
bien_raw <- BIEN_occurrence_species(species_names,all.taxonomy=T,native.status=T,
                                  natives.only=F,observation.type=T,collection.info=T,political.boundaries=T,
                                  cultivated=T)
nrow(bien_raw) #2326
# write file
write.csv(bien_raw, "bien_raw.csv")

bien_raw <- read.csv("bien_raw.csv")

################################################################################
# F) USDA Forest Service, Forest Inventory and Analysis (FIA) download
################################################################################

# First, download raw data
# Go to https://apps.fs.usda.gov/fia/datamart/CSV/datamart_csv.html
# Either download the "TREE" file (e.g., "AL_TREE.csv") for each state
#   (works well if you only need a few) or scroll to the bottom of the
#   page and download "TREE.csv", which gives data for all states combined
#   (9.73 GB); you need lots of memory to do it with just the one "TREE" file
# Place all the tree files in an "fia_read_in" folder in your working
#   directory
# While you're on the FIA data download webpage, scroll to the bottom of
#   the page and download the "PLOT.csv" file (data for all states combined)
#   and place in your working directory

# read in FIA species codes
setwd("G:/My Drive/Hoban_Lab_Docs/Projects/Butternut_JUCI/DataFiles/Occurrence/fia_read_in")
fia_code <- 0601
# join taxa list to FIA species codes
species_code <- 0601

# I have to read in each state file separately and pull info for our target
#   taxa then remove the file before loading the next state because memory
#   gets exhaused otherwise



# function to extract target species data from each state CSV
extract_tree_data <- function(file_name){
  # read in tree data, which lists all species and the plots in which they were
  #   found; larger ones will take time to read in
  state_df <- read.csv(file_name)
  # cycle through vector of target species codes and extract those rows from
  #   the state CSV
  for (sp in 1:length(fia_code)){
    target_sp <- state_df[which(state_df$SPCD==fia_code[[sp]]),]
    fia_raw <- rbind(fia_raw,target_sp)
  }
  # remove state file to make space for reading in next one
  rm(state_df)
  # take a look at how much data were pulled
  print(dim(fia_raw))
  return(fia_raw)
}

# make a new data frame to gather data for target taxa
fia_raw <- data.frame()
# create list of state files
file_list <- list.files(path = "G:/My Drive/Hoban_Lab_Docs/Projects/Butternut_JUCI/DataFiles/Occurrence/fia_read_in/FIA_states",
                        pattern = ".csv",full.names = T)
# loop through states and pull data using the function defined above
fia_outputs <- lapply(file_list, extract_tree_data)

# stack state-by-state data extracted to create one dataframe
fia_raw <- data.frame()
for(file in seq_along(fia_outputs)){
  fia_raw  <- rbind(fia_raw , fia_outputs[[file]])
}; nrow(fia_raw) #3055
# write file
write.csv(fia_raw,"fia_raw.csv")


### DESCRIPTION:
# This script provides instructions and code chunks for downloading wild
#   occurrence points from:
# GLOBAL DATABASES (though all likely have U.S. bias?)
# Global Biodiversity Information Facility (GBIF)
# Integrated Digitized Biocollections (iDigBio)
# U.S. Herbarium Consortia (SERNEC, SEINet, etc.)
# Botanical Information and Ecology Network (BIEN)
# NATIONAL DATABASES
# Forest Inventory and Analysis (FIA) Program of the USDA Forest Service

### INPUT:
# target_taxa_with_syn.csv (list of target taxa)
# columns:
# 1. "taxon_name" (genus, species, infra rank, and infra name, all
#    separated by one space each; hybrid symbol should be " x ", rather
#    than "_" or "???", and go between genus and species)
# 2. (optional) "taxon_name_acc" (accepted taxon name you have chosen)
# 3+ (optional) other data you want to keep with taxa info

### OUTPUTS:
# gbif_raw.csv
# idigbio_raw.csv
# herbaria_raw.csv
# bien_raw.csv
# fia_raw.csv

#################
### LIBRARIES ###
#################

library(plyr)
library(tidyverse) #ggplot2,dplyr,tidyr,readr,purrr,tibble,stringr,forcats
library(data.table)
library(batchtools)
library(textclean)


setwd("G:/My Drive/Hoban_Lab_Docs/Projects/Butternut_JUCI/DataFiles/Occurrence")

bien_raw <- read.csv("bien_raw.csv")
idigbio_raw <- read.csv("idigbio_raw.csv")
fia_raw <- read.csv("G:/My Drive/Hoban_Lab_Docs/Projects/Butternut_JUCI/DataFiles/Occurrence/fia_read_in/fia_raw.csv")
gbif_raw <- read.csv("G:/My Drive/Hoban_Lab_Docs/Projects/Butternut_JUCI/DataFiles/Occurrence/gbif_read_in/gbif_raw.csv")
sernec_raw <- read.csv("G:/My Drive/Hoban_Lab_Docs/Projects/Butternut_JUCI/DataFiles/Occurrence/sernec_read_in/occurrences.csv")

#################
### FUNCTIONS ###
#################

# searches for data frame columns with only NAs and removes them
remove.empty.col <- function(df){
  remove <- vector(mode = "character")
  for(i in 1:ncol(df)){
    if(sum(is.na(df[,i])) == nrow(df)){
      remove <- c(remove,names(df)[i])
      print(names(df)[i])
    }
  }
  if(length(remove)>0){
    df <-  df[,-which(names(df) %in% remove)]
  }
  return(df)
}

fia_raw <- remove.empty.col(fia_raw)

write.csv(fia_raw, "fia_raw.csv")

# calculates percent of each data frame column that is not NA
percent.filled <- function(df){
  for(i in 1:ncol(df)){
    print(paste(names(df)[i],": ",
                round((nrow(df)-sum(is.na(df[,i])))/nrow(df),3)*100,"%",sep=""))
  }
}

################################################################################
# 1. Read in target taxa list and raw occurrence data
################################################################################

setwd("./../..")
setwd("G:/My Drive/Hoban_Lab_Docs/Projects/Butternut_JUCI/DataFiles/Occurrence")

gbif_raw <- read.csv("gbif_read_in/gbif_raw.csv",header=T,
                     na.strings=c("","NA"),stringsAsFactors=F)
gbif_raw <- remove.empty.col(gbif_raw) #; percent.filled(gbif_raw)
idigbio_raw <- read.csv("idigbio_raw.csv",header=T,
                        na.strings=c("","NA"),stringsAsFactors=F)
idigbio_raw <- remove.empty.col(idigbio_raw) #; percent.filled(idigbio_raw)
sernec_raw <- read.csv("sernec_read_in/occurrences.csv",header=T,
                       na.strings=c("","NA"),stringsAsFactors=F)
sernec_raw <- remove.empty.col(sernec_raw) #; percent.filled(sernec_raw)
bien_raw <- read.csv("bien_raw.csv",header=T,
                     na.strings=c("","NA"),stringsAsFactors=F)
bien_raw <- remove.empty.col(bien_raw) #; percent.filled(bien_raw)
fia_raw <- read.csv("fia_read_in/fia_raw.csv",header=T,
                    na.strings=c("","NA"),stringsAsFactors=F)
fia_raw <- remove.empty.col(fia_raw) #; percent.filled(fia_raw)

################################################################################
# 2. Subset and standardize column names
################################################################################

# GBIF
gbif_raw$taxon_name <- "Juglans cinerea"
gbif_raw$dataset <- "GBIF"
# create taxon_name column

gbif_df <- data.frame(cbind(as.character(gbif_raw$identifier), gbif_raw$taxon_name, gbif_raw$dataset, gbif_raw$decimalLongitude, gbif_raw$decimalLatitude))
colnames(gbif_df) <- c("ID", "TaxonName","DataSet","Longitude","Latitude")
  
# iDigBio

# split eventdate column to just get year

idigbio_raw$taxon_name <- "Juglans cinerea"
idigbio_raw$dataset <- "iDigBio"
idigbio_df <- data.frame(cbind(idigbio_raw$catalognumber, idigbio_raw$taxon_name, idigbio_raw$dataset, idigbio_raw$geopoint.lon,idigbio_raw$geopoint.lat))

colnames(idigbio_df) <- c("ID", "TaxonName","DataSet","Longitude","Latitude")

# US_Herbaria
# create taxon_name column
sernec_raw$taxon_name <- "Juglans cinerea"
sernec_raw$dataset <- "US_Herbarium"

sernec_df <- data.frame(cbind(sernec_raw$catalogNumber, sernec_raw$taxon_name, sernec_raw$dataset, sernec_raw$decimalLongitude, sernec_raw$decimalLatitude))
colnames(sernec_df) <- c("ID", "TaxonName","DataSet","Longitude","Latitude")

# BIEN
bien_raw$taxon_name <- "Juglans cinerea"
bien_raw$dataset <- "BIEN"
bien_df <- data.frame(cbind(bien_raw$catalog_number, bien_raw$taxon_name, bien_raw$dataset, bien_raw$longitude, bien_raw$latitude))
colnames(bien_df) <- c("ID", "TaxonName","DataSet","Longitude","Latitude")

# FIA

# read in supplemental FIA tables:
#  - list of species tracked and their codes
#  - state and county codes and names
#  - plot level data (has lat-long)
setwd("G:/My Drive/Hoban_Lab_Docs/Projects/Butternut_JUCI/DataFiles/Occurrence/fia_read_in/FIA_states")

fia_codes <- 0601
county_codes <- read.csv("US_state_county_FIPS_codes.csv", header = T,
                         na.strings=c("","NA"), colClasses="character")
#county_codes <- read.csv("US_state_county_FIPS_codes.csv")
fia_plots <- read.csv("PLOT.csv")
# remove unnecessary columns from plot data
fia_plots <- fia_plots[,c("INVYR","STATECD","UNITCD","COUNTYCD","PLOT",
                          "LAT","LON")]
# join FIA data to supplemental tables

fia_raw2 <- join(fia_raw,fia_plots)
fia_raw2$taxon_name <- "Juglans cinerea"
fia_raw2$dataset <- "FIA"

fia_df <- data.frame(cbind(fia_raw2$CN, fia_raw2$taxon_name, fia_raw2$dataset, fia_raw2$LON, fia_raw2$LAT))

colnames(fia_df) <- c("ID","TaxonName","DataSet","Longitude","Latitude")

##

all_occurrences_butternut <- rbind(bien_df, fia_df, gbif_df, idigbio_df, sernec_df)

##cleanup 

butternut_complete_occurrence <- all_occurrences_butternut[complete.cases(all_occurrences_butternut), ]
write.csv(butternut_complete_occurrence, "butternut_occurrence.csv")

##plot

newmap <- getMap(resolution = "low")
plot(newmap)
plot(newmap, xlim = c(max(as.character(butternut_complete_occurrence$Longitude)),min(as.character(butternut_complete_occurrence$Longitude)), ylim = c(max(as.character(butternut_complete_occurrence$Latitude)),min(as.character(butternut_complete_occurrence$Latitude)))
points(as.character(butternut_complete_occurrence$Longitude), as.character(butternut_complete_occurrence$Latitude), col = "dodgerblue", cex = 0.8, pch = 16)


##setwd
setwd("G:/My Drive/Hoban_Lab_Docs/Projects/Butternut_JUCI/DataFiles/Occurrence")

##plot occurrence records

butternut_oc_matrix <- as.matrix(butternut_complete_occurrence[,4:5])


colnames(butternut_oc_matrix) <- c("Longitude", "Latitude")

worldmap <- getMap(resolution = "low")
setwd("G:/My Drive/Hoban_Lab_Docs/Projects/Butternut_JUCI/DataFiles/Occurrence")
pdf("butternut_occurrence.pdf", width = 10, height = 8)
plot(newmap, xlim = c(-109.55523,9.566805), ylim = c(-37,73.966667))

points(as.character(butternut_complete_occurrence[,4]),as.character(butternut_complete_occurrence[,5]), pch =  16, col = "dodgerblue")
dev.off()

##
plot(newmap, xlim = c(min(butternut_latlon$Longitude),max(butternut_latlon$Longitude)), ylim = c(min(butternut_latlon$Latitude),max(butternut_latlon$Latitude)))




setwd("G:/My Drive/Hoban_Lab_Docs/Projects/Butternut_JUCI/DataFiles/44Populations")

##then calculate mean latitude for each population
butternut_latlon <- read.csv("Butternut_44pops_latlon.csv") 
butternut_mean_lon <- matrix()
butternut_mean_lat <- matrix()

##identifying mean latitudes for each population

for(pop in butternut_44pop_names){
  
  butternut_mean_lon[pop] <- mean(butternut_latlon[butternut_latlon$Pop == pop,][,3])
  
  
}

for(pop in butternut_44pop_names){
  
  butternut_mean_lat[pop] <- mean(butternut_latlon[butternut_latlon$Pop == pop,][,4])
  
  
}

##convert to matrix
butternut_mean_lon <- matrix(butternut_mean_lon)
butternut_mean_lat <- matrix(butternut_mean_lat)

##document cleanup
butternut_mean_lon <- butternut_mean_lon[-1]
butternut_mean_lat <- butternut_mean_lat[-1]

#combine into one document for mean long and lat for each pop
butternut_coords <- matrix(ncol = 2, nrow = 44)
butternut_coords[,1] <- butternut_mean_lon
butternut_coords[,2] <- butternut_mean_lat
rownames(butternut_coords) <- butternut_44pop_names
colnames(butternut_coords) <- c("Mean Lon", "Mean Lat")

##plot

pdf("butternut_occurrence_zoom.pdf", width = 10, height = 8)
butternutmap <- getMap(resolution = "low")
plot(butternutmap, xlim = c(min(butternut_coords[,1]), max(butternut_coords[,1])), ylim = c(min(butternut_coords[,2]),max(butternut_coords[,2])))
points(butternut_oc_matrix[,1],butternut_oc_matrix[,2],pch = 16, col = "dodgerblue")
dev.off()

