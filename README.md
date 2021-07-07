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

## Uppset plot
![upset_locations](https://user-images.githubusercontent.com/25608100/124742211-c610cc00-df1c-11eb-82a1-28b4f61600cc.png)

It's a graph to visualize intersections of multiple sets compared to the traditional approaches, i.e. the Venn Diagram.

#### Prepare the data from pyloseq object
```{.r}
#Delete the singeltons.
physeq_single <- filter_taxa(physeq_object, function (x) {sum(x > 0) > 1}, prune=TRUE)

#Extract the transpose ASV's table
physeq_asvtable<-t(otu_table(physeq_single))

#Extract the metadata
physeq_metadata <- data.frame(phyloseq::sample_data(physeq_single), 
                       check.names = FALSE
		       )

#Create the ASV's list from each group that we want show.
ASVlist_physeq_groupA = colnames(physeq_asvtable[physeq_metadata$Variable == "A", 
			apply(physeq_asvtable[physeq_metadata$Variable == "A",], MARGIN=2, function(x) any(x >0))])

#Create the data for the graph
data_upset = list(A=ASVlist_physeq_groupA,B=ASVlist_physeq_groupB)

#Calculate all combinations
combinations_upset = make_comb_mat(data_upset)
```
#### Create the graph

Variables to change:
* 120 <- minimum of interactions to show.
* ylim <- te minimum and max of the total asv by sample.
```{.r}
#Create the graph
upset_graph<-UpSet[comb_size(combinations_upset) >= 5],comb_order = rev(order(comb_size(combinations_upset[comb_size(
	    combinations_upset) >= 120]))),right_annotation = upset_right_annotation(combinations_upset,
	    ylim = c(0, 6000)),width = unit(ifelse(1,20), "cm") )

##Paramters of the graph
ss_upset = set_size(combinations_upset[comb_size(combinations_upset) >= 5])
cs_upset = comb_size(combinations_upset[comb_size(combinations_upset) >= 5])
od_upset = column_order(upset_graph)

##Make graph
png("graph/upset_graph.png",height=1200,width=2600,res=150)
ht = draw(upset_graph)
decorate_annotation("Intersection\nsize", {
    grid.text(cs_upset[od_upset], x = seq_along(cs_upset), y = unit(cs_upset[od_upset], "native") + unit(4, "pt"), 
        default.units = "native", just = c("center", "bottom"), 
        gp = gpar(fontsize = 8, col = "#404040"), rot = 0, )
})
dev.off()
```

