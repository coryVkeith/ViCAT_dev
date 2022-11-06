### Author: Cory Von Keith
### Email: coryvonkeith@gmail.com
### OrcID: https://orcid.org/0000-0002-7111-3497

###This code is to be used with data outputs from ViCAT and will result in the generation of 
###relative abundance bar plots, relative abundance heat maps, and ordination maps to compare
###between samples.


#Install Packages if needed.

install.packages('installr')
require(installr)
install.packages('compositions')
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("phyloseq")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")

#Load Packages

library(tidyverse)
library(phyloseq); packageVersion("phyloseq")
library("ggplot2")
library("plyr")
library("dplyr")
library('compositions')
library(DESeq2);packageVersion("DESeq2")

#___________________________________________________________________________________________|
#                                                                                           |
#Before you move to this step, manually alter the metadata.csv file to include relevant info|
#for your project.                                                                          |
#___________________________________________________________________________________________|



#___________________________________________________________________________________________|
#                                                                                           |
###Normalization of data is an important step for accurate comparison of communities between|
###samples. This normalization uses DeSeq2 for median of ratios normalization. You may wish |
###to have a different normalization step, so change this step if needed.                   |
#                                                                                           |
#___________________________________________________________________________________________|

mydata <- read_csv("ViCAT_otutable.csv", col_names = TRUE)

norm_data <- read_csv("ViCAT_otutable.csv", col_names = TRUE)
norm_data
norm_data <- as.matrix(norm_data[,-1])

meta <- read_csv("metadata.csv", col_names = TRUE)
meta <- meta %>% 
  column_to_rownames(var="PlantID")
meta

#Checks if files have same columns and rows
all(colnames(norm_data) %in% rownames(meta))
all(colnames(norm_data) == rownames(meta))

dds <- DESeqDataSetFromMatrix(countData = norm_data, colData = meta, design = ~ Genotype)
dds
View(counts(dds))

dds <- estimateSizeFactors(dds, type = 'poscounts')
#dds <- estimateSizeFactors(dds)
sizeFactors(dds)

###Normalize
normalized_counts <- counts(dds, normalized=TRUE)
write.table(normalized_counts, file="ViCAT_normalized_counts.tsv", sep="\t", quote=F, col.names=NA)

###Replacing the Taxa assignment
normalized_counts <- read_tsv("ViCAT_normalized_counts.tsv", col_names = TRUE)
normalized_counts

### Takes the column from mydata (Accession numbers)
refid <- mydata[1]
refid

###Replaces normalized data with the Accession numbers and assigns the column name with RefID
normalized_counts[1] <- refid
colnames(normalized_counts) [1] <- "RefID"
normalized_counts

mydata <- normalized_counts                                                                         
mydata                                                                                            

#___________________________________________________________________________________________________|
#                                                                                                   |
#                         This next block creates a Phyloseq Object.                                |
#                                                                                                   |
#___________________________________________________________________________________________________|

### create otu matrix for read counts
otumat <-  mydata %>%
  column_to_rownames(var="RefID")
otumat <- as.matrix(otumat)
colnames(otumat) <- paste0("Sample", 1:ncol(otumat))


### create taxa matrix with species names
taxmat <- read_tsv("sim_mock_taxmat.tsv") %>%
  column_to_rownames(var="RefID")
taxmat$Species <- as.character(taxmat$Species)
taxmat <- as.matrix(taxmat)


### create meta data matrix
metamat <- read_csv("metadata.csv", col_names = TRUE)
metamat <- metamat %>% 
  column_to_rownames(var="SampleID")

####Variables for physeq object
META = sample_data(metamat, errorIfNULL = TRUE)
OTU = otu_table(otumat, taxa_are_rows = TRUE)
TAX = tax_table(taxmat)

###Create phyloseq object
physeq = phyloseq(OTU, TAX, META)

#___________________________________________________________________________________________________|
#                                                                                                   |
#                         This next block creates a generates figures for analysis.                 |
#                                                                                                   |
#___________________________________________________________________________________________________|


###Transforms data to compositional data, and plots Relative Abundance Bar Plot
ra= transform_sample_counts(physeq, function(x) x/sum(x))
barplot <- plot_bar(ra, "PlantID", fill = "Species")
png(file="Barplot.png", width=1000, height=1000)
plot(barplot)
dev.off()
###Creates heatmap
heatmap <- plot_heatmap(ra, "NMDS", "jaccard", "PlantID")
png(file="Heatmap_jaccard.png", width=500, height=500)
plot(heatmap)
dev.off()

                            
###### This code block will plot all distance methods

###Calls distance methods available for use.
dist_methods <- unlist(distanceMethodList)
print(dist_methods)

###Sets the distance methods to use 
dist_methods <- dist_methods[-(1:3)]
dist_methods["designdist"]

###List of plots that the following for loop will write to.
plist <- vector("list", length(dist_methods))
names(plist) = dist_methods

###Iterates through the distance methods and creates MDS plots for each.
for( i in dist_methods ){
  # Calculate distance matrix
  iDist <- phyloseq::distance(physeq, method=i)
  # Calculate ordination
  iMDS  <- ordinate(physeq, "MDS", distance=iDist)
  ## Make plot
  # Don't carry over previous plot (if error, p will be blank)
  p <- NULL
  # Create plot, store as temp variable, p
  p <- plot_ordination(physeq, iMDS)
  # Add title to each plot
  p <- p + ggtitle(paste("MDS using distance method ", i, sep=""))
  # Save the graphic to file.
  plist[[i]] = p
}
### Plots all of the generated plots
df = ldply(plist, function(x) x$data)
names(df)[1] <- "distance"
p = ggplot(df, aes(Axis.1, Axis.2, color=Genotype))
p = p + geom_point(size=3, alpha=0.5)
p = p + facet_wrap(~distance, scales="free")
p = p + ggtitle("MDS using various distance metrics for dataset")
p
png(file="MDS_dist_metrics.png", width=1000, height=1000)
plot(p)
dev.off()                            
                            
                            
                            
                            
#___________________________________________________________________________________________|
#                                                                                           |
#The following code will generate NMDS and PCoA ordination plots using the desired distance |
# metrics as determined by evaluating the data with the MDS using all distances from the    |
# code above. Change the other variables for desired plot style and title. Other options    |
# can be changed depending on desired aesthetics.                                           |
#___________________________________________________________________________________________|                            
                            
### Change this to the desired distance metric.                           
dist = "jaccard"
### Change this to have the colors of the points in the ordination match the desired column in the metadata file.
COLOR = "Genotype"
### Change this to change the plot title.
TITLE = "Sim_Mock"
                            
###NMDS 
dist = "jaccard"
ordi = ordinate(physeq, method="NMDS", distance=dist)
capture.output(ordi, file = "NMDS_stress.txt", append = FALSE)
NMDS = plot_ordination(physeq, ordi, "samples", color=COLOR)
NMDS = NMDS + ggtitle(TITLE) + geom_polygon(aes(fill=COLOR))
NMDS
png(file="NMDS.png", width=500, height=500)
plot(NMDS)
dev.off()

###PCoA 
physeq.ord <- ordinate(physeq, "PCoA", distance=dist, binary = FALSE)
PCoA = plot_ordination(physeq, physeq.ord, type="samples", color=COLOR) 
PCoA = PCoA + geom_polygon(aes(fill=COLOR)) + geom_point(size=5) + ggtitle(TITLE)
PCoA
png(file="PCoA.png", width=500, height=500)
plot(PCoA)
dev.off()



