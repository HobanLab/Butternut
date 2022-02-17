<b><p><h1 style="color:red;font-size:20px;"> Project Description</b></p></h1>
This project file contains R Scripts created for a project on <i>Juglans cinerea</i>, also known as butternut. We are interested in determining how butternut's range has shifted in response to modern and past climate changes. Hoban et al. (2010) determined that butternut's genetic patterns were most inluenced by post-glacial migration northward, so this study added individuals sampled by the Jeanne Romero-Severson and Martin William's labs to identify if these results are consistent when adding more comprehensive sampling in this species' northern range. We also wanted to more specifically examine past sources of data to identify the range shifts of this taxon. Therefore the main focuses on this project are listed here: 

![Alt text](Images/worflow_github.jpg?raw=true "description of the multi-faceted approach to this project, using genetic analyses, hindcast species distribution models, and fossil pollen records to determine how butternut's past range shifts influenced current genetic diversity.") 

The code written for performing two of the steps are in this Github; code for downloading, projecting, and creating pollen and glacier maps were created by Alissa Brown (https://github.com/alissab/juglans). This repository is split into multiple folders: Archive, Images, SDMs, and Genetic_Analyses.

<b><p><h1 style="color:red;font-size:20px;">Folder Descriptions</b></p></h1>

<b>Archive:</b> This is a folder for code that was used in initial steps of this analysis but did not end up in the final manuscript. Most of this code was performed on more individuals than ended up in the final results (were removed due to isolation or genetic relatedness) or for loci that differed based on scoring year. 

<b>Images:</b> This folder contains images generated for organizing the Github and population maps. 

<b>SDMs:</b> This folder is for the code to run "species distribution models," specifically with boosted regression tree models as described in Elith et al. (2008). The code published here is based on code from Peter Breslin (ASU) and Fabio Suzart de Albuquerque. Part of this project was interested in identifying butternut's ecological preferences using species distribution modeling and then modeling habitat requirements into the past. This folder contains the code for performing all of these analyses. Here is a conceptual diagram of the steps and files names: 

![Alt text](Images/SDM_flowchat.jpg?raw=true "Flowchart for species distribution models, further described in text below")  

<p>The order these steps were performed in is indicated with the number in their name. All R Scripts used for SDM analysis have "SDM" written next the number and are briefly summarized below:</p>
<ol start="1">
<li>First, occurrence records must be downloaded from herbarium databases (in this case, Global Biodiversity Information Facility (GBIF), Integrated Digitized Biocollections (iDigBio), South East Regional Network of Expertise and Collections (SERNEC), Botanical Information and Ecology Network (BIEN), and the national network of forest survey plots managed by the Forest Inventory and Analysis Program (FIA) of the USDA Forest Service) using code that was based off code created by Emily Beckman (https://github.com/esbeckman/IMLS_Beckman). Then, these occurrences were used to create a range polygon that buffered out 100 km from the range edges (which was eventually used for genetic diversity analyses).</li>
<li> Occurrence records were reduced for spatial autocorrelation by selecting one point for every point that is within 1 km of others.</li>
<li>Pseudo-absence points were generated in equal number to presence points. </li>
<li>Following this, 19 Worldclim variables were extracted at every point. Variables that were most correlated with presence and least correlated with one another (correlation coefficient < 0.5) were selected to model species distribution. This was done by generating a dissimilarity matrix and choosing the variable that had the highest correlation with presence in each cluster of highly correlated climate variables.</li>
<li>The final boosted regression tree model was run with five variables and from these predictors a model of suitable habitat was made and then hindcast into past climate scenarios using eight time periods representing notable periods in post-LGM climatic history, available from the Paleoclim database (Brown et al., 2018): 130 (last interglacial); 22 (LGM); 17-14.7 (Heinrich-Stadial); 14.7-12.9 (Bolling-Allerod); 12.9-11.7 (Younger Dryas); 11.7-8.326 (early Holocene); 8.326-4.2 (mid-Holocene); and 4.2-0.3 (late Holocene)  ka YBP (thousand years before present). </li>
</ol>

<b> Files within SDMs folder</b> 
<ul><li>InputFiles</li></ul>
<ul><ul><li>Description: All input files utilized for generating butternut's species distribution model, habitat suitability maps, and projecting the habitat suitability maps into the past.</li></ul></ul>
<ul><ul><li>bio_2-5m_bil</li></ul></ul>
 <ul><ul><ul><li>Description: Folder containing all 19 current (averages over 1970 - 2013) bioclimatic layers used to generate butternut's species distribution model, downloaded at 2.5 m resolution. </li></ul></ul></ul>
 
 <ul><ul><li>bound_p</li></ul></ul>
 <ul><ul><ul><li>Description: Folder containing the shapefile used to limit the area of extent to make butternut range maps. </li></ul></ul></ul>
 
 <ul><ul><li>occurrence_records</li></ul></ul>
 <ul><ul><ul><li>Description: The folder containing all of the raw records used to create the full occurrence record file for this species - titled, butternut_complete_occurrence. </li></ul></ul></ul>

 <ul><ul><li>Paleo_Files</li></ul></ul>
 <ul><ul><ul><li>Description: Within this file is the values of all 19 bioclimatic variables during the last 130,000 years. These variables were used to project butternut's distribution model into the past and are separated into 8 separate time points. These bioclimatic values were all downloaded at 2.5 m resolution from Paleoclim.</li></ul></ul></ul>
 
<ul><ul><li>butternut_abs: CSV of pseudo-absence points used in generating the boosted regression trees that is used to predict butternut's stuitable habitat.</li></ul></ul>
<ul><ul><li>butternut_pa: CSV of all presence and pseuod-absence points used in the BRT model generating butternut's species distribution, with an additional column coded indicating what is a presence point (1) and what is an absence point (0). </li></ul></ul>
<ul><ul><li>butternut_var: CSV of all presence and pseudo-absence points (indicated by a 1 and 0, respectively) with the values of all 19 bioclimatic variables at the location of the occurrence record extracted to each point. </li></ul></ul>
<ul><ul><li>elevation_extent: TIF file of North American elevation limited to the extent of the analysis.</li></ul></ul>
<ul><ul><li>extent_project: TIF file of North American elevation limited to the extent of this analysis project to Albers Equal Area Conic projection.</li></ul></ul>
<ul><ul><li>occurrence_noauto_noproj: Occurrence records used in this analysis, cleaned for spatial autocorrelation and not projected.</li></ul></ul>
 
<ul><li>OutputFiles</li></ul>
<ul><ul><li>Description: All of the result files generated when creating the species distribution model of butternut.</ul></ul></li>
<ul><ul><li>PaleoClim_HSMs</li></ul></ul>
<ul><ul><ul><li>Description: PDF and images of the habitat suitability maps for butternut during the last 130,000 to present. </li></ul></ul></ul>
<ul><ul><li>biseral_cor_matrix: CSV of the biserial correlation coefficients between presence points and each bioclimatic variable.</ul></ul></li>
<ul><ul><li>contribution_bp_allvar: PDF of selected bioclim variables used to generate butternut's SDM and their percent contribution to the final SDM. </ul></ul></li>
<ul><ul><li>contribution_bp_allvar: PDF and image of the bioclim variables dissimilarity neighbor-joining tree. Variables are joined by least dissimilar to most dissimilar. Threshold for non-significance autocorrelation is marked with a red line at 0.5 ranked dissimilarity.</ul></ul></li>
<ul><ul><li>extent_pointsnoauto: PDF of butternut's habitat suitability map with occurrence records plotted on it. </ul></ul></li>
<ul><ul><li>hsm: PDF of butternut's habitat suitability map. </ul></ul></li>
<ul><ul><li>hsm_a: PDF of butternut's habitat suitability map with pseudo-absence points plotted on top of it. </ul></ul></li>
<ul><ul><li>hsm_pa: PDF of butternut's habitat suitability map with presence and pseudo-absence points plotted on top of it. </ul></ul></li>
 <ul><ul><li>hsm_worldclim_only_raster: Raster of butternut's habitat suitability map. </ul></ul></li>
  <ul><ul><li>worldclim_stack: RData file, TIF file, and GRD file of all selected worldclim variables stacked. </ul></ul></li>
 
<b>Genetic Analyses: </b> The genetic diversity portion of this project is contained in the genetic_analyses folder, which contains the R Scripts to run genetic diversity and structure analyses, along with the regressions between genetic diversity and geographic location. We used genetic data from the publication Hoban et al. (2010) and newer sampling efforts on butternut from 2011 - 2015. These individuals were collected by Jeanne Romero-Severson, Sean Hoban (https://github.com/smhoban), and Martin Williams over the course of near ten years with a major sampling effort closer to 2009 and then followed up by another round of sampling 2012 - 2015. The initial individuals that were collected were genotyped by Sean Hoban and then subsequent individuals were genotyped in the Romero-Severson lab at Notre Dame non-consequetively. The order these analyses were performed in is indicated with a preceeding number. The script labeled "comparison_barplot" was code used to determine if there were scoring differences between researcher which led to some re-binning analyses to ensure consistency of allele scoring when researchers differed. Then, the code for removing individuals based on missing data and relatedness was designed so PCoA and structure could be run. The individuals were also plotted on a map following removal for missing data. Also, mean latitude, longitude, allelic richness and heterozygosity were calculated for all populations. Finally, a loop was written to compare mean latitude and distance to range edge of each population to genetic diversity.

![Alt text](Images/gendiv_flowchart.jpg?raw=true "Flowchart for the storage of data files. R Scripts to run genetic analyses are contained within this folder and they are all run on the data files stored in the data_files folder. The results are then stored in the genetic_analyses_results folder.") 

<b>Folders within Genetic_Analyses</b> 

data_files: Data files generated throughout the process of some of the analyses are stored here. Before_reorg contains the 3 population genind used for the comparison barplot code, after_reorg contains all of the genind and data frames generated following the re-binning process. geographic_files contains many of the geographic files used in analyses, including the range shapefile with the 100 km buffer (butternut_buffer) and many of the mea population longitude and latitude data frames. 

genetic_analyses_results: All of the outputs from genetic analyses that were used in the final publication, including figures and data frames, are stored here. 

<b>Files within Genetic_Analyses</b>
<ul><li>data_files</li></ul>
<ul><ul><li>after_reorg</li></ul></ul>
<ul><ul><ul><li>Description: Following geographic reorganization into populations and rebinning analysis, these are the files used in the final genetic analyses included in the manuscript. </li></ul></ul></ul>
<ul><ul><ul><li> butternut_24pop.gen: Genepop file of 1721 butternut individuals following geographic reorganization into 24 populations. </li></ul></ul></ul>
<ul><ul><ul><li>butternut_24pop_nomd: Genepop, Genalex, CSV, and Arlequin files of genotypes of butternut individuals used for main genetic analyses (cleaned for clones and missing data). When lonlat is included, it includes coordinates and when "relate" is included, its used for relatedness analysis. </li></ul></ul></ul>
<ul><ul><ul><li>butternut_24pop_relate_red: CSV, genepop, and Arlequin files of genotypes of butternut individuals, reduced for individuals with missing data and too high of relatedness, leaving 993 individuals. </li></ul></ul></ul> 
<ul><ul><ul><li>butternut_26pop_nomd: Genepop and Arlequin files of genotypes of butternut individuals including Quebec population individuals. CSV file titled "lonlat" includes coordinate information for all individuals. </li></ul></ul></ul>
<ul><ul><ul><li>butternut_26pop_relate_red: Genepop, Arlequin, and CSV files of genotypes of butternut individuals reduced by relatedness but with Quebec butternut individuals. </li></ul></ul></ul>

<ul><ul><li>before_reorg</li></ul></ul>
<ul><ul><ul><li>Description: Genepop files from initial analyses on butternut individuals; rebinning analysis and geographic analyses were performed on these individuals to then yield the data files in the "after_reorg" file. </li></ul></ul></ul> 
<ul><ul><ul><li>butternut_3pops: Genepop file of all initial 1761 butternut individuals, sorted by the person who scored data files. This document was used in the initial rebinning analysis to determine if there were differences in scoring with different individual scorers. </li></ul></ul></ul> 
<ul><ul><ul><li>butternut_44pops: Genepop file of all initial 1,761 which was used for initial genetic and geographic analyses.</li></ul></ul></ul> 

<ul><ul><li>geographic_files</li></ul></ul>
<ul><ul><ul><li>Description: Files used in generating maps and all geographic analyses for butternut. </li></ul></ul></ul> 
<ul><ul><li>butternut_buffer: Shapefile of butternut's total range, generated from occurrence records and buffered out by 100 km and cleaned to have smooth edges.</li></ul></ul>
<ul><ul><ul><li>buffernut_buffer_map: Map of mean population coordinates plotted over the butternut range buffer.</li></ul></ul></ul>
<ul><ul><ul><li>butternut_coord_df: CSV of the mean coordinates of the main 24 populations used in this analysis with raw names and color-coded by population.</li></ul></ul></ul>
<ul><ul><ul><li>butternut_dist_edge_df: CSV of the distance to range edge for each of the 24 populations.</li></ul></ul></ul>
<ul><ul><ul><li>max_min_lonlat: CSV of the minimum and maximum coordinates of the butternut range.</li></ul></ul></ul>
 
<ul><li>genetic_analyses_results</li></ul>
<ul><ul><li>Clustering_Analyses</li></li></ul>
<ul><ul><ul><li>STR</li></ul></ul></ul>
<ul><ul><ul><ul><li>Butternut_Structure_24pops: Folder containing results of structure runs on 993 individuals.</li></ul></ul></ul></ul>
<ul><ul><ul><ul><li>Butternut_Structure_26pops: Folder containing results of structure runs on 1005 individuals./li></ul></ul></ul></ul>
<ul><ul><ul><li>Butternut_24pop_PCoA: PDF, PNG, and Inkscape file of a PCoA of 993 butternut individuals.</li></ul></ul></ul>
<ul><ul><ul><li>Butternut_24pop_Structure: Structure diagram of 993 individuals of the best supported K value for structure. </li></ul></ul></ul>
<ul><ul><ul><li>Butternut_24pop_Structure_alldeltaK_Supplement.PNG: Image file of several delta K values for 993 butternut individuals, presented in the supplement.</li></ul></ul></ul>
<ul><ul><ul><li>Butternut_26pop_PCoA_Structure: PDF, PNG, and Inkscape file of a PCoA of the 1005 butternut individuals including Quebec individuals.</li></ul></ul></ul>
 
<ul><ul><li>Diversity_Analyses</li></li></ul>
<ul><ul><ul><li>butternut_24pop_allrich_rp_df.csv: R2 and p-values for regressions between each 24 butternut populations' geographic statistics (mean population latitude; distance to range edge) and allelic richness. </li></ul></ul></ul>
<ul><ul><ul><li>butternut_24pop_gendiv_stat_df.csv: Genetic diversity summary stats of 24 populations of butternut populations.</li></ul></ul></ul>
<ul><ul><ul><li>butternut_24pop_hexp_rp_df.csv: R2 and p-values for regressions between each 24 butternut populations' between geographic statistics (mean population latitude; distance to range edge) and expected heterozygosity. </li></ul></ul></ul>
<ul><ul><ul><li>butternut_24pop_hwe.csv: Hardy-Weinberg Equilibrium deviation by population 24 populations.</li></ul></ul></ul>
<ul><ul><ul><li>butternut_24pop_ld_loci.csv: Linkage disequilibrium between each of the 11 loci. </li></ul></ul></ul>
<ul><ul><ul><li>butternut_26pop_allrich_rp_df</li></ul></ul></ul>
<ul><ul><ul><li>butternut_26pop_hexp_rp_df</li></ul></ul></ul>

<ul><ul><li>Rebinning_Analyses</li></li></ul>
<ul><ul><ul><li>PCoA</li></ul></ul></ul>

<b><p><h1 style="color:red;font-size:20px;">R Information</b></p></h1>

R version 4.0.5 (2021-03-31)
Platform: x86_64-w64-mingw32/x64 (64-bit)
Running under: Windows 10 x64 (build 19041)

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
[1] poppr_2.9.2      Demerelate_0.9-3 adegenet_2.1.3   ade4_1.7-17    diveRsity_1.9.90    rgdal_1.5-23 sp_1.4-5 pegas_1.0-1        ape_5.5            PopGenReport_3.0.4 hierfstat_0.5-7 poppr_2.9.2     

<b><p><h1 style="color:red;font-size:20px;">References</b></p></h1>
Brown, J. L., Hill, D. J., Dolan, A. M., Carnaval, A. C., & Haywood, A. M. (2018). PaleoClim, high spatial resolution paleoclimate surfaces for global land areas. Scientific data, 5(1), 1-9

Elith, J., Leathwick, J. R., & Hastie, T. (2008). A working guide to boosted regression trees. Journal of Animal Ecology, 77(4), 802-813.

Fick, S. E., & Hijmans, R. J. (2017). WorldClim 2: new 1‐km spatial resolution climate surfaces for global land areas. International journal of climatology, 37(12), 4302-4315.

Hoban, S. M., Borkowski, D. S., Brosi, S. L., McCLEARY, T. S., Thompson, L. M., McLachlan, J. S., ... Romero-Severson, J. (2010). Range‐wide distribution of genetic diversity in the North American tree Juglans cinerea: A product of range shifts, not ecological marginality or recent population decline. Molecular Ecology, 19(22), 4876-4891.
