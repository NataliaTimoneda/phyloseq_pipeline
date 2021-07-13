This pipeline contains a pipeline for a metagenomic analysis using phyloseq package from R.


## Libraries used in al pipeline
<details><summary>All the libraries</summary>
<p>
	
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
</p>
</details>

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
If you need to delete samples by the metadata. 

```{.r}
physeq_v2 = subset_samples(physeq_object, column_name != "variable")
```

## Alpha & Beta diversity

In this example we want to compare the alpha diversity between the variable "Depth".
In this example we compare the betadiversity and mark each sample with 3 variables.

<img src="https://user-images.githubusercontent.com/25608100/125424544-a06d6c55-cfcc-462d-9c7a-789eff2d7beb.png" width="600" />

<details><summary> The code</summary>
<p>

```{.r}
## Ordre of the Variable
newSTorder = c("0", "5", "10","20", "40","60", "75")

## Make the pairs comparisons, with depth variable
compare_depth <- combn(levels(sample_data(physeq_object)[,"Depth"]$Depth)[-9],2,
                    simplify=FALSE)

## Alpha diversity by Depth and Sampling
# If you want to add several samples in each Depth category change the variable in the parameter "x".

plot<-plot_richness(physeq_object, color = "Depth", x="Depth" ) + geom_boxplot()
plot$data$Depth <- factor(plot$data$Depth, levels=newSTorder)

## Create the plot
png("graph/alfa_diversity_depth.png",height=1200,width=2000,res=150)
plot
dev.off()
```
</p>
</details>

#### Beta diversity

The input is a phyloseq object with data normalize and delete the samples with no ASV's. (This happens when you select a specific taxonomy)

<img src="https://user-images.githubusercontent.com/25608100/125425991-770d7280-fc86-4de0-94b0-af6561b5bde3.png" width="600" />

<details><summary> The code</summary>
<p>

You can change the color and the shape of each point, and add a label.

```{.r}
bray_phyloseq = ordinate(phyloseq_normalize, "NMDS", "bray", weighted = TRUE)
stressplot(bray_phyloseq)

p.bray = plot_ordination(phyloseq_normalize, bray_phyloseq, color = "Variable1", shape= "Variable2",
 	 title = "NMDS of Bray-Curtis distance") + geom_text(mapping = aes(label = Varible3 ), size = 3, vjust = 1.5)
```
</p>
</details>

## Uppset plot

It's a graph to visualize intersections of multiple sets compared to the traditional approaches, i.e. the Venn Diagram.

Complete inforamation: https://jokergoo.github.io/ComplexHeatmap-reference/book/upset-plot.html

![upset_locations](https://user-images.githubusercontent.com/25608100/124742211-c610cc00-df1c-11eb-82a1-28b4f61600cc.png)

<details><summary> The code</summary>
<p>

#### Prepare the data from pyloseq object
```{.r}
#Delete the singeltons. Delete the ASV's only have 1 sequence.

physeq_single <- filter_taxa(physeq_object, function (x) {sum(x > 0) > 1}, prune=TRUE)

#Extract the transpose ASV's table
physeq_asvtable<-t(otu_table(physeq_single))

#Extract the metadata
physeq_metadata <- data.frame(phyloseq::sample_data(physeq_single), 
                       check.names = FALSE
		       )

#Create the ASV's list from each group that we want show. One for each group:
ASVlist_physeq_groupA = colnames(physeq_asvtable[physeq_metadata$Variable == "A", 
			apply(physeq_asvtable[physeq_metadata$Variable == "A",], MARGIN=2, function(x) any(x >0))])

#To create automatically the previous step:
#Create a vector with all the group of the variable
list<-levels(metadata_elu$Varible)

#Make one variable for each group
for (group in list) {
    name<-paste("ASVlist_phyloseq_", group, sep="")
    assign(name, colnames(physeq_asvtable[physeq_metadata$Variable == group, apply(physeq_asvtable[physeq_metadata$Variable == group,], MARGIN=2, function(x) any(x >0))]))
}

#Create the data for the graph
data_upset = list(A=ASVlist_physeq_groupA,B=ASVlist_physeq_groupB)

#To create automatically the previous step:
data_upset<-""
for (group in list) {
    name<-paste("ASVlist_phyloseq_", group, sep="")
    to_add= paste(group,"=",name, sep="")
    data_upset<-c(data_upset,to_add )
}

#Calculate all combinations
combinations_upset = make_comb_mat(data_upset)
```
#### Create the graph

Variables to change:
* 120 <- minimum of interactions to show.
* ylim <- te minimum and max of the total asv by sample.

```{.r}
#Create the graph
upset_graph<-UpSet[comb_size(combinations_upset) >= 120],comb_order = rev(order(comb_size(combinations_upset[comb_size(
	    combinations_upset) >= 120]))),right_annotation = upset_right_annotation(combinations_upset,
	    ylim = c(0, 6000)),width = unit(ifelse(1,20), "cm") )

##Paramters of the graph
ss_upset = set_size(combinations_upset[comb_size(combinations_upset) >= 120])
cs_upset = comb_size(combinations_upset[comb_size(combinations_upset) >= 120])
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

</p>
</details>

