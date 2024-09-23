# Cancer Genomics Analyses
# By PARADA, CA

# Clean up environment
rm(list = ls())

# Set your working directory
setwd("C:\\Users\\cdsp\\Desktop\\Parada\\")

# Install readbulk package
# install.packages("readbulk")  

# Load readbulk
library(readbulk)


################################################################################
################################################################################
# In your local computer,
## create a directory named raw_data_TN
### add all tumor-normal SAMPLE_combined_indels_PASS.hg38_multianno.txt and SAMPLE_combined_snvs_PASS.hg38_multianno.txt in there

## create a directory named "raw_data_TO"
## add all 
# I do recommend using the "normal" pipeline we have running in the lab
# The high_confidence pipeline may be filtering out too many variants at this point; therefore,
# it requires additional optimization
################################################################################


################################################################################
# Tumor Normal Outputs
################################################################################
# Read and combine files. The files will be merged into one dataframe.
merge_TN <- read_bulk(directory = "C:\\Users\\cdsp\\Desktop\\Parada\\Chordoma_Project\\WES_Outputs\\tumor-normal\\raw_data_TN",
                   fun=read.delim, 
                   stringsAsFactors=FALSE)

# Let's look at the column names
colnames(merge_TN)

# Note that, in the merge_TN file, the column "File" contains the entire sample name
merge_TN$File

# Let's keep only the sample ID in the column "File" for clarity
merge_TN$sample <- merge_TN$File # create a new column named "sample" and duplicate column "File" 
# check
colnames(merge_TN) 
merge_TN$sample

# Split column "sample" into "sampleID" and "fileextension"
newCols <- c('SET', 'sampleID', 'combined', 'Vtype', 'PASS', 'multianno')
# Check
merge_TN <- separate(merge_TN, col = 'sample', into = newCols, sep = '_')

# Check
unique(merge_TN$SET)
unique(merge_TN$sampleID)
unique(merge_TN$Vtype)

# combine the prefix 'SET' to the sampleID 
library(dplyr)
# combine the prefix SET to sampleID to follow the same sampleID prefix seen in the merge_TN dataset
merge_TN$sampleID <- paste(merge_TN$SET, merge_TN$sampleID, sep = "_")
merge_TN$sampleID

# Delete cols that we are not going to use
colnames(merge_TN)
merge_TN[, c(199:200, 202)] <- NULL
merge_TN[,201:202] <- NULL
colnames(merge_TN)

# Add column mode indicating TN pipeline
merge_TN$mode <- 'TN'
unique(merge_TN$mode)

# Save file in CSV
write.csv(merge_TN, 'merge_TN.csv', row.names = F)

################################################################################
# Tumor Only outputs
################################################################################
merge_TO <- read_bulk(directory = "C:\\Users\\cdsp\\Desktop\\Parada\\Chordoma_Project\\WES_Outputs\\tumor-only\\raw_data_TO",
                      fun=read.delim, 
                      stringsAsFactors=FALSE)

# Let's look at the column names
colnames(merge_TO)

# Note that, iIn the merge file, the column "File" contains the entire sample name
merge_TO$File

# Let's keep only the sample ID in the column "File" for clarity
merge_TO$sample <- merge_TO$File # create a new column named "sample" and duplicate column "File" 
# check
colnames(merge_TO) 
merge_TO$sample
# Split column "sample" into "sampleID" and "fileextension"
newCols <- c('SET', 'sampleID', 'mode', 'mutect2', 'somatic', 'filtered', 'Vtype', 'normalized', 'PASS', 'multianno')
# Check
merge_TO <- separate(merge_TO, col = 'sample', into = newCols, sep = '_')
# Check
unique(merge_TO$sampleID)
unique(merge_TO$Vtype)
unique(merge_TO$SET)
# combine the prefix 'SET' to the sampleID 
library(dplyr)
# combine the prefix SET to sampleID to follow the same sampleID prefix seen in the merge_TN dataset
merge_TO$sampleID <- paste(merge_TO$SET, merge_TO$sampleID, sep = "_")
merge_TO$sampleID

# Delete cols that we are not going to use
colnames(merge_TO)
merge_TO[, c(198:199, 202:204, 206:208)] <- NULL
colnames(merge_TO)

# Save file in CSV
write.csv(merge_TO, 'merge_TO.csv', row.names = F)

################################################################################
# Creating a maf object
################################################################################
# Install maftools if needed
#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

#BiocManager::install("maftools")

# Load maftools 
library(maftools)

# Convert the annotated VCF into a maf object using annovarToMaf()
TN <- annovarToMaf('C:\\Users\\cdsp\\Desktop\\Parada\\merge_TN.csv',
             tsbCol = "sampleID",
             table = "refGene",
             ens2hugo = TRUE,
             basename = NULL,
             sep = ",",
             MAFobj = F,
             sampleAnno = NULL,
             refBuild = "hg38")

TO <- annovarToMaf('C:\\Users\\cdsp\\Desktop\\Parada\\merge_TO.csv',
                   tsbCol = "sampleID",
                   table = "refGene",
                   ens2hugo = TRUE,
                   basename = NULL,
                   sep = ",",
                   MAFobj = F,
                   sampleAnno = NULL,
                   refBuild = "hg38")

multimerge <- rbind(TO, TN, fill = T)

# Save multimerge file
write.csv(multimerge, 'multimerge.csv', row.names = F)

################################################################################
# Filtering maf
################################################################################
# keep only exonic mutations
multimergeFiltered <- multimerge[which(multimerge$Func.refGene == 'exonic'), ]
unique(multimergeFiltered$Func.refGene) # check if there are only exonic variants

# Exclude synonymous SNV and unknown variants
multimergeFiltered <- multimergeFiltered[which(multimergeFiltered$ExonicFunc.refGene != 'synonymous SNV'), ]
multimergeFiltered <- multimergeFiltered[which(multimergeFiltered$ExonicFunc.refGene != 'unknown'), ]
unique(multimergeFiltered$ExonicFunc.refGene)# check if synonymous SNV and unknown variants are deleted

#### Split the values of column Otherinfo13 
# The column `multimergeFiltered$Otherinfo13` contains the values of genotype (GT), Allele Depth (AD), Allele Frequency (AF), Depth (DP), etc all combined in this single column, separated by `:`
# Let's organize these values storing them in single columns.
multimergeFiltered$Otherinfo12 # variables
multimergeFiltered$Otherinfo13 # values
newCols <- c('GT','AD','AF','DP','F1R2','F2R1','FAD','PGT','PID','PS','SB')
multimergeFiltered <- separate(multimergeFiltered, col = 'Otherinfo13', into = newCols, sep = ':')
multimergeFiltered$AF # check if now there is a separated column for AF and DP
multimergeFiltered$DP
multimergeFiltered$AD

#### identifying and labeling tumor suppressor/oncogenes
# Load our gene list of known tumor suppressor/oncogenes
OGTS <- read.csv('CP-OGTS-list.csv', header = T)
# Note that our gene list have 2 columns: Hugo_Symbol and is.driver == 'Yes'
colnames(OGTS)

# Load dplyr package to use the join function in r
library(dplyr)

# use the join function to joicombine the list of suppressor genes or oncogenes (potential drivers) with your dataset
multimergeFiltered <- left_join(multimergeFiltered, OGTS, by = 'Hugo_Symbol')
# Save formatted dataset
write.csv(multimergeFiltered, "multimergeFilteredOGTS.csv", row.names = F)


################################################################################
# Maftools analysis
################################################################################
# Source: https://bioconductor.org/packages/release/bioc/vignettes/maftools/inst/doc/maftools.html

library(maftools)
# Summarizing maf
# Make a vector containing all values for variant classification
# This is going to be a list of variant classifications to be considered as non-synonymous.
vc_nonSyn <- multimergeFiltered$Variant_Classification
unique(vc_nonSyn) # check the list created

# Create a dataframe only with known altered cancer drivers (is.driver == 'Yes')
mutatedDrivers <- multimergeFiltered[which(multimergeFiltered$is.driver == 'Yes'), ]

# Create a maf object using the function read.maf()
maf <- read.maf(multimergeFiltered, vc_nonSyn = vc_nonSyn)
mafDrivers <- read.maf(mutatedDrivers, vc_nonSyn = vc_nonSyn)
# Note that: vc_nonSyn parameter is NULL by default.
# We are manually providing a list of variant classifications to be considered as non-synonymous. 
# Rest will be considered as silent variants. 
# Default uses Variant Classifications with High/Moderate variant consequences.


# clinical information, if any, can be uploaded and included in this step. 
# clinicalData <- read.csv('path/to/clinical/data/clinicalData.csv', header = T, sep = "\,") 
# maf <- read.maf(multimergeFiltered, vc_nonSyn = vc_nonSyn, clinicalData = clinicalData)

#Shows sample summary.
getSampleSummary(maf)
#Shows gene summary.
getGeneSummary(maf)
#shows clinical data associated with samples
getClinicalData(maf)
#Shows all fields in MAF
getFields(maf)
# Write all maf summaries to your working dir under the basename of your choice
write.mafSummary(maf = maf, basename = 'maf')

## Visualization of maf summaries
#We can use the plotmafSummary() function to plot the summary of the maf file, 
#which displays number of variants in each sample as a stacked barplot and 
#variant types as a boxplot summarized by Variant_Classification.
plotmafSummary(maf = maf, 
               rmOutlier = TRUE, 
               addStat = 'median', 
               dashboard = TRUE, # If FALSE plots simple summary instead of dashboard style.
               titvRaw = FALSE,
               fs = 1, # base size for text
               textSize = 0.5, # font size if showBarcodes is TRUE.
               showBarcodes = TRUE, # include sample names in the top bar plot.
               top = 20) # include n top genes. 

plotmafSummary(maf = mafDrivers, 
               rmOutlier = TRUE, 
               addStat = 'median', 
               dashboard = TRUE, # If FALSE plots simple summary instead of dashboard style.
               titvRaw = FALSE,
               fs = 1, # base size for text
               textSize = 0.5, # font size if showBarcodes is TRUE.
               showBarcodes = TRUE, # include sample names in the top bar plot.
               top = 20) # include n top genes. 

# From now on, run the analysis with both maf and mafDrivers
# Visualize results side-by-side and compare them

#### Visualization of top 25 altered cancer driver genes
# Gene list
genes <- maf@data[which(maf@data$is.driver == 'Yes'), ]
genes <- genes$Hugo_Symbol

# use the function mafbarplot() to visualize the top 25 altered genes 
mafbarplot(maf,
           n = 25,
           color = NULL,
           fontSize = 0.7,
           includeCN = FALSE,
           legendfontSize = 0.7,
           borderCol = "#34495e",
           showPct = TRUE)

## Visualization of maf file by oncoplot (waterfall plot)
# Better representation of maf file can be shown as oncoplots, also known as waterfall plots.
# Visualize Top50 genes
oncoplot(maf = maf, 
         showTumorSampleBarcodes = T,
         SampleNamefontSize =1.3,
         fontSize = 0.4,
         n = 25)

# Visualize drivers
oncoplot(maf = mafDrivers, 
         showTumorSampleBarcodes = T,
         SampleNamefontSize =1.3,
         fontSize = 0.4,
         top = 25)


# Visualize genes of interest
geneOfInterest <- c("ATRX","ATM", "DOCK1", "GNAS", "JAG2", "LAMA2", "MAP3K21", "MAP3K6", "PIK3C2G")

oncoplot(maf = maf, 
         showTumorSampleBarcodes = T,
         SampleNamefontSize =0.8,
         fontSize = 0.6,
         genes = geneOfInterest)
# Note that variants annotated as Multi_Hit are those genes which are mutated more than once in the same sample.

# You can also usp the oncotrip() to simply the visualization of geneOfInterest
oncostrip(maf = maf, 
          showTumorSampleBarcodes = T,
          SampleNamefontSize =0.8,
          fontSize = 0.6,
          genes = geneOfInterest)


# Rainfall plots

# Cancer genomes, especially solid tumors are characterized by genomic loci with localized hyper-mutations
# Such hyper mutated genomic regions can be visualized by plotting inter variant distance on a linear genomic scale. 
# These plots are called rainfall plots and we can draw such plots using rainfallPlot. 
# If the parameter detectChangePoints = TRUE, rainfall plot also highlights regions where 
# potential changes in inter-event distances are located.

# Plot only shows the results for the most mutated sample
rainfallPlot(maf =  maf,
             detectChangePoints = TRUE, # If TRUE it writes kataegis results as a tab delimted file.
             pointSize = 0.4, 
             ref.build = "hg38",
             fontSize = 0.7)
?rainfallPlot()
## Visualizing Variant Allele Frequency (VAF)
# The function `plotVaf()` plots Variant Allele Frequencies as a boxplot
# This quickly helps to estimate clonal status of top mutated genes 
# Clonal genes usually have mean allele frequency around ~50% assuming pure sample.
plotVaf(maf = maf, vafCol = 'AF')
plotVaf(maf = mafDrivers, vafCol = 'AF')

## Somatic Interactions
# Mutually exclusive or co-occurring set of genes can be detected using the `somaticInteractions()` function

# Visualize exclusive/co-occurrence events on top 10 mutated genes. 
somaticInteractions(maf = maf,
                    top = 10, 
                    pvalue = c(0.05, 0.1),
                    genes = NULL, # You can assess interactions among genes of interest only (i.e genes = geneOfInterest)
                    fontSize = 0.8,
                    countsFontSize = 0.8,
                    sigSymbolsSize = 2,
                    sigSymbolsFontSize = 0.9)

somaticInteractions(maf = mafDrivers,
                    top = 10, 
                    pvalue = c(0.05, 0.1),
                    genes = NULL, # You can assess interactions among genes of interest only (i.e genes = geneOfInterest)
                    fontSize = 0.8,
                    countsFontSize = 0.8,
                    sigSymbolsSize = 2,
                    sigSymbolsFontSize = 0.9)

## Predicting cancer driver genes based on positional clustering
# The function `oncodrive()` identifies disease drivers from the maf.
# oncodrive is based on the oncodriveCLUST Python algorithm. 
# Concept is based on the fact that most of the variants in cancer causing genes are enriched at few specific loci (aka hot-spots). 
# This method takes advantage of the hotspots to identify cancer drivers. 
# If you use this function, please cite [Tamborero2013](https://academic.oup.com/bioinformatics/article/29/18/2238/240376)
dev.off()
maf.sig = oncodrive(maf = maf, AACol = 'aaChange', minMut = 3, pvalMethod = 'zscore')

# The function `plotOncodrive()` shows a scatter plot with size of the points proportional to the number of clusters found in the gene. 
plotOncodrive(res = maf.sig, fdrCutOff = 0.05, useFraction = F, labelSize = 0.7, bubbleSize = 4)

## Visualize pfam Domains
# The function pfamDomains() adds pfam domain information to the amino acid changes of your variants
# pfamDomain also summarizes amino acid changes according to the domains that are affected. 
#This serves the purpose of knowing what domain in given cancer cohort, is most frequently affected. 
dev.off()
laml.pfam = pfamDomains(maf = maf, AACol = 'aaChange', top = 10)
?pfamDomains()

## Visualize Drug-Gene Interactions
# The function drugInteractions() checks for drug-gene interactions and gene druggability 
# based on the Drug gene Interactions database
# You can assess the drug interactions and druggability of the top altered genes of your dataset
dgi = drugInteractions(maf = maf, fontSize = 1.2, plotType = 'bar', top = 50, drugs = T)
dgiDrivers = drugInteractions(maf = mafDrivers, fontSize = 1.2, plotType = 'bar', top = 50, drugs = T)
# or only for the list of genes of interest
dgiInterest = drugInteractions(maf = maf, fontSize = 1.2, plotType = 'bar', genes = geneOfInterest, drugs = T)
# save table
write.csv(dgi, 'dgi.csv', row.names = F)
write.csv(dgiDrivers, 'dgiDrivers.csv', row.names = F)
write.csv(dgiInterest, 'dgiInterest.csv', row.names = F)

#### Lolliplot visulization of mutations and hotspots of interest 
# the `lollipopPlot()` function requires amino acid changes information in the maf file.
# This information should be displayed in the aaChange column of your MAF file.

#lollipop plot for DNMT3A, which is one of the most frequent mutated gene in Leukemia.
geneOfInterest <- c("ATRX","ATM", "DOCK1", "GNAS", "JAG2", "LAMA2", "MAP3K21", "MAP3K6", "PIK3C2G")
dev.off()

# plotting one gene
lollipopPlot(maf = maf, 
             gene = 'GNAS', 
             AACol = 'aaChange',
             labelPos = 'all',
             labPosSize = 1.2,
             domainLabelSize = 0.8)
# plotting multiple genes
lollipopPlot(maf = maf, 
             gene = c('GNAS', 'JAG2'), 
             AACol = 'aaChange',
             labelPos = 'all',
             labPosSize = 0.7,
             domainLabelSize = 0.9)

#### Another way to plot lolliplots
# Another way to plot lolliplots is using the bioconductor package `trackviewer`

# Install trackviewer
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("trackViewer")

# Load trackviewer
library(trackViewer)

# For example, let's plot the ABCC2 mutations on the positions S281N and  A985A
# First, let's set up the gene name and the mutation(s) found on this gene 
protein <- GRanges("ABCC2", 
                   IRanges(start=c(281, 535), # mutation position in the protein without the 1 letter amino acid abbreviation
                           width = 1, 
                           names = c("S281N", "A985A"))) # mutation position in the protein with the 1 letter amino acid abbreviation 
# Note that the arguments for the parameter `start` and `names` can be found in the aaChange column of your MAF file

# Set up gene name and length 
features_protein <- GRanges("ABCC2", # gene name
                            ranges = IRanges(start =0, end = 1545)) # protein length. You can look it up on uniprot.org 

# Now you can plot the lolliplot based on the information you entered for the protein and features_protein
lolliplot(protein, features_protein)

# Color mutation lollipops
protein$color <- c("red", "red") 

# Color protein sequence
features_protein$fill <- "purple"

# Visualize the customized plot
lolliplot(protein, features_protein)

#### General protein domains can be drawn with the function plotProtein()
# The function `plotProtein()` automatically plots the domains of the protein of interest.
# All you need is is the gene name and respective unique identifier in the NCBI RefSeq database
plotProtein(gene = "TP53", refSeqID = "NM_000546")

sessionInfo()




