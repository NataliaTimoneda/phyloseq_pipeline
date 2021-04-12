---
title: "Phyloseq pipeline"
\author: Natàlia Timoneda
geometry: margin=1.5cm
output: pdf_document
---
This pipeline contains a pipeline for a metagenomic analysis using phyloseq package from R.


## Libraries used in al pipeline

```{.r}
library(ggforce)
library(fantaxtic)
library(yingtools2)
library(dplyr)
library(ape)
library(reshape)
library(phyloseq)
library(ggpubr)
library(data.table)
library(vegan)
library(pairwiseAdonis)
library(ggplot2)
library(reshape2)
library(grid)
library(gridExtra)
library(ranacapa)
library(devtools)
source_url("https://raw.githubusercontent.com/obigriffith/biostar-tutorials/master/Heatmaps/heatmap.3.R")
library("gplots")
library("devtools")
library("DESeq2")
library("ampvis")
library("ampvis2")
library(metagMisc)
library(VennDiagram)
library("corrplot")
```

## Load the data
```{.r}
OTUtable_file <- read.table("data/asv_table.tbl", row=1, header=TRUE,sep="\t")
OTUtable <- otu_table(OTUtable_file, taxa_are_rows=TRUE)
file_metadata <- read.csv("data/metadata.csv", header=TRUE, row.names=1, sep=",")
metadata = sample_data(data.frame(file_metadata))
file_taxa <- read.table("data/taxonomy.tbl", header=TRUE, row.names=1) 
taxa_tmp <- as.matrix(file_taxa)
TAXtable = tax_table(taxa_tmp)

physeq_object  = merge_phyloseq(OTUtable,metadata,TAXtable)
```
If ypu need to delete samples by the metadata. 

```{.r}
physeq_v2 = subset_samples(physeq, column_name != "variable")
```


## Alpha & Beta diversity
```{.r}
## Ordre of the depth variable
newSTorder = c("0", "5", "10","20", "40","60", "75")

## Make the pairs comparisons, with depth variable
compare_depth <- combn(levels(sample_data(physeq_v2)[,"Depth"]$Depth)[-9],2,
                    simplify=FALSE)

## Alpha diversity by Depth and Sampling
p<-plot_richness(physeq_v2, color = "Depth", x="Sampling" ) + geom_boxplot()
p$data$Depth <- factor(p$data$Depth, levels=newSTorder)

## Create the plot
png("graph/alfa_diversity_all_depth_sample.png",height=1200,width=2000,res=150)
p
dev.off()
```