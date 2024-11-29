####    ****************load libraries************    ##### 
library(tidyverse)
library(DESeq2)
library(GEOquery)
library(tidyverse)
library(gridExtra)
library(WGCNA)
library(CorLevelPlot)
library(ggplot2)
library(dplyr)
library(tibble)
library(org.Hs.eg.db)
library(clusterProfiler)
library(RColorBrewer)
library(dendextend)
library(pheatmap)
library(igraph) 
library(RColorBrewer) 

allowWGCNAThreads()  


# **********load Counts file***********   #
counts_data <- read.csv(file = "Counts_Adult.csv", header = TRUE, sep = ",")
counts_data  <- as.data.frame(counts_data)
counts_data <- column_to_rownames(counts_data , var ="GeneID")

data <- counts_data


# ****metadata****
sample_info <- read.csv(file = "sampleInfo_Adult.csv", header = TRUE, sep = ",")

#--------------------------------------------------------------------------------------------------------------------------------


# *****create dds******
dds <- DESeqDataSetFromMatrix(countData = data,
                              colData = sample_info,
                              design = ~ 1)


dds75<- dds[rowSums(counts(dds) >= 25) >= 15,]  # Filtering
nrow(dds75) #14433 genes 

# ***performing variance stabilization **** #
dds_norm <- vst(dds75)

#get normalized values
norm.counts <- assay(dds_norm) %>%
  t()

if ("THAP9" %in% colnames(norm.counts)) {
  print("Column exists!")
} else {
  print("Column does not exist.")
}

#--------------------------------------------------------------------------------------------------------
# PCA plot to remove outliers 
#--------------------------------------------------------------------------------------------------------
pcaData <- plotPCA(dds_norm, intgroup = c("cell_type"), returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))


# Plot PCA *************
ggplot(pcaData, aes(x = PC1, y = PC2, color = cell_type)) +
  geom_point(size = 3) +
  geom_text(aes(label = rownames(pcaData)), vjust = -0, hjust = -0) +  
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed() +
  theme_bw()

#------------------------------------------------------------------------------------------------------------------
# Dendrogram with colored samples to remove the outliers
#-------------------------------------------------------------------------------------------------------------------
norm.counts <- norm.counts %>%
  t()

# ***hierarchical clustering**  #
htree_norm <- hclust(dist(t(norm.counts)), method = "average")


cell_types <- unique(sample_info$cell_type)
colors <- brewer.pal(length(cell_types), "Dark2")  
names(colors) <- cell_types

sample_colors_norm <- colors[sample_info$cell_type]
dend_norm <- as.dendrogram(htree_norm) # Converting to dendrogram


ordered_sample_colors_norm <- sample_colors_norm[order.dendrogram(dend_norm)]

labels_colors(dend_norm) <- ordered_sample_colors_norm

###   **********Save the plot**************   ###
pdf("dendrogram_plot_norm_adult_new.pdf", width = 12, height = 8)
par(mar = c(5, 4, 4, 2) + 0.1)
plot(dend_norm, main = "Hierarchical Clustering of Adult Samples")
dev.off()

#---------------------------------------------------------------------------------------------------------------
# choosing a set of soft-thresholding powers for network construction
#---------------------------------------------------------------------------------------------------------------------

power <- c(c(1:10), seq(from = 12, to = 50, by = 2))

norm.counts <- t(norm.counts)

# to call the network topology analysis function 
sft <- pickSoftThreshold(norm.counts,
                         powerVector = power,
                         networkType = "signed",
                         verbose = 5)

sft.data <- sft$fitIndices

####        ***visualization to pick power***      ####
a1 <- ggplot(sft.data, aes(Power, SFT.R.sq, label = Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  geom_hline(yintercept = 0.8, color = 'red') +
  labs(x = 'Power', y = 'Scale free topology model fit, signed R^2') +
  theme_classic()

a2 <- ggplot(sft.data, aes(Power, mean.k., label = Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  labs(x = 'Power', y = 'Mean Connectivity') +
  theme_classic()

grid.arrange(a1, a2, nrow = 2)


#------------------------------------------------------------------------------------------------------------------------
# blockwise net module function 
#--------------------------------------------------------------------------------------------------------------------------------

# converting matrix to numeric, in order to use it for blockwise function 
norm.counts[] <- sapply(norm.counts, as.numeric)
soft_power <- 8
adjacency = adjacency(norm.counts, power = soft_power)
TOM = TOMsimilarity(adjacency)  

temp_cor <- cor
cor <- WGCNA::cor

# memory estimate w.r.t blocksize
bwnet <- blockwiseModules(norm.counts,
                          maxBlockSize = 9000,
                          minModuleSize = 30,
                          TOMType = "signed",
                          power = soft_power,
                          mergeCutHeight = 0.25,
                          numericLabels = FALSE,
                          randomSeed = 1234,
                          verbose = 3)

cor <- temp_cor

#---------------------------------------------------------------------------------------------------------------
# Module eigengenes and plotting of dendrogram with merged and unmerged colors
#--------------------------------------------------------------------------------------------------------
module_eigengenes <- bwnet$MEs
head(module_eigengenes)

table(bwnet$colors) # To get number of genes for each module

# Plot the dendrogram and the module colors before and after merging 
PDC1 <- plotDendroAndColors(bwnet$dendrograms[[1]], 
                            cbind(bwnet$unmergedColors[bwnet$blockGenes[[1]]], bwnet$colors[bwnet$blockGenes[[1]]]),
                            c("unmerged", "merged"),
                            dendroLabels = FALSE,
                            addGuide = TRUE,
                            hang = 0.03,
                            guideHang = 0.05)


PDC2 <- plotDendroAndColors(bwnet$dendrograms[[2]], 
                            cbind(bwnet$unmergedColors[bwnet$blockGenes[[2]]], bwnet$colors[bwnet$blockGenes[[2]]]),
                            c("unmerged", "merged"),
                            dendroLabels = FALSE,
                            addGuide = TRUE,
                            hang = 0.03,
                            guideHang = 0.05)

#--------------------------------------------------------------------------------------------------
# table of which gene is in which module using blockwise net module function 
#--------------------------------------------------------------------------------------------------

genes_modules <- bwnet$colors 

genes_modules <- as.data.frame(genes_modules)

genes_modules <- genes_modules %>%
  rownames_to_column(var= "GeneID")


#-------------------------------------------------------------------------------------------------
# Defining levels and binarizing categorical variable  
#-------------------------------------------------------------------------------------------------

rownames(sample_info) <- sample_info$sample_id 
sample_info$cell_type <- factor(sample_info$cell_type, levels = c("MO", "OPC"))

# Creating binary columns using model.matrix
binary_columns <- model.matrix(~ cell_type - 1, data = sample_info)
traits <- binary_columns

# Defining numbers of genes and samples
nSamples <- nrow(norm.counts)
nGenes <- ncol(norm.counts)


#------------------------------------------------------------------------------------------------------
# module trait correlation
#------------------------------------------------------------------------------------------------------

# calculating correlation between module eigengenes and traits 
module.trait.corr <- cor(module_eigengenes, traits, use = 'p')
module.trait.corr.pvals <- corPvalueStudent(module.trait.corr, nSamples)

##  **visualizing using of heat map**   ## 
heatmap.data <- merge(module_eigengenes, traits, by = 'row.names')
head(heatmap.data)

heatmap.data <- heatmap.data %>% 
  column_to_rownames(var = 'Row.names')

textMatrix <- paste(signif(module.trait.corr, 2), "\n(", signif(module.trait.corr.pvals, 1), ")", sep = "")
dim(textMatrix) <- dim(module.trait.corr)

pheatmap(module.trait.corr,
         cluster_rows = FALSE, 
         cluster_cols = FALSE,
         display_numbers = textMatrix,
         fontsize_number = 7,
         fontsize_row = 13,
         fontsize_col = 13,
         angle_col = 0,
         main = "Module Trait Correlation Heatmap",
         color = colorRampPalette(c("blue1", "white", "red"))(50))

#----------------------------------------------------------------------------------------------------------------
#Module membership and gene significance with specific trait calculation 
#----------------------------------------------------------------------------------------------------------------

modNames = substring(names(module_eigengenes), 3)

#module membership ****
module.membership.measure <- cor(module_eigengenes, norm.counts, use = 'p')
module.membership.measure.pvals <- corPvalueStudent(module.membership.measure, nSamples)

module.membership.measure <- as.data.frame(module.membership.measure)
module.membership.measure.pvals <- as.data.frame(module.membership.measure.pvals)

# Calculate the gene significance and associated p-values ******
traits <- as.data.frame(traits)

gene.signf.corr.MO <- cor(norm.counts, traits$`cell_typeMO`, use = 'p')
gene.signf.corr.pvals.MO <- corPvalueStudent(gene.signf.corr.MO, nSamples)

gene.signf.corr.OPC <- cor(norm.counts, traits$`cell_typeOPC`, use = 'p')
gene.signf.corr.pvals.OPC <- corPvalueStudent(gene.signf.corr.OPC, nSamples)

#-----------------------------------------------------------------------------------------------------------------
# to plot the module membership(MM) vs gene trait significance (GS) of specific module 
#------------------------------------------------------------------------------------------------------------------
modulenames <- substring(names(module_eigengenes),3)
modulecolors <- bwnet$colors
modulecolors <- as.data.frame(modulecolors)

module = "turquoise"
column = match(module, modulenames)
moduleGenes = modulecolors==module

module.membership.measure <- module.membership.measure %>%
  t()

# Extract the x and y coordinates for THAP9 for MO
x_coord <- abs(module.membership.measure["THAP9", column])
y_coord <- abs(gene.signf.corr.MO["THAP9", 1])

sizeGrWindow(7,7)
par(mfrow = c(1,1))
verboseScatterplot(abs(module.membership.measure[moduleGenes, column]),
                   (abs(gene.signf.corr.MO[moduleGenes,1])),
                   xlab = paste("Module Membership of Turquoise Module"),
                   ylab = "Gene significance for Mature Oligodendrocyte",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2 , col = module)
points(x_coord, y_coord, pch = 10, col = "red", cex = 1.5)
text(x_coord, y_coord + 0.03, labels = "THAP9", col = "red", cex = 1.2, font = 2)

# Extract the x and y coordinates for THAP9 for OPC
y_coord_OPC <- abs(gene.signf.corr.OPC["THAP9", 1])


sizeGrWindow(7,7)
par(mfrow = c(1,1))
verboseScatterplot(abs(module.membership.measure[moduleGenes, column]),
                   (abs(gene.signf.corr.OPC[moduleGenes,1])),
                   xlab = paste("Module Membership of Turquoise Module"),
                   ylab = "Gene significance for Oligodendrocyte Progenitor cells",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2 , col = module)
points(x_coord, y_coord_OPC, pch = 10, col = "red", cex = 1.5)
text(x_coord, y_coord_OPC + 0.03, labels = "THAP9", col = "red", cex = 1.2, font = 2)

#---------------------------------------------------------------------------------------------------------------------
# genes with high MM vs GS (Hub Genes) from the turquoise module cut-off 0.8 
#--------------------------------------------------------------------------------------------------------------------
modulecolors <- bwnet$colors
module = "turquoise"

turquoise_genes <- names(modulecolors)[modulecolors == module]


MM_turquoise_0.8 <- module.membership.measure[turquoise_genes, column][abs(module.membership.measure[turquoise_genes, column]) >= 0.8] # 0.8 is threshold
MM_turquoise_0.8 <- as.data.frame(MM_turquoise_0.8)
GS_turquoise_0.8_MO <- gene.signf.corr.MO[turquoise_genes, 1, drop = FALSE][abs(gene.signf.corr.MO[turquoise_genes, 1]) >= 0.8, , drop = FALSE]
GS_turquoise_0.8_MO <- as.data.frame(GS_turquoise_0.8_MO)


common_hubgenes_MO <- intersect(rownames(MM_turquoise_0.8), rownames(GS_turquoise_0.8_MO))

# Creating a data frame for the genes in the turquoise module
hubgenes_turquoise_0.8_MO <- data.frame(
  Gene = common_hubgenes_MO,
  Module_Membership = module.membership.measure[common_hubgenes_MO, column],
  Gene_Significance = gene.signf.corr.MO[common_hubgenes_MO, 1])


GS_turquoise_0.8_OPC <- gene.signf.corr.OPC[turquoise_genes, 1, drop = FALSE][abs(gene.signf.corr.OPC[turquoise_genes, 1]) >= 0.8, , drop = FALSE]
common_hubgenes_OPC <- intersect(rownames(MM_turquoise_0.8), rownames(GS_turquoise_0.8_OPC))

hubgenes_turquoise_0.8_OPC <- data.frame(
  Gene = common_hubgenes_OPC,
  Module_Membership = module.membership.measure[common_hubgenes_OPC, column],
  Gene_Significance = gene.signf.corr.OPC[common_hubgenes_OPC, 1])


all(rownames(GS_turquoise_0.8_OPC) %in% rownames(GS_turquoise_0.8_MO))
all(rownames(GS_turquoise_0.8_OPC) == colnames(GS_turquoise_0.8_MO))


#-------------------------------------------------------------------------------------
# Edge and node files for network visualization 
#--------------------------------------------------------------------------------------

# export the gene list of new module *****
for (i in 1:length(bwnet$MEs)){
  modules = c(substring(names(bwnet$MEs)[i], 3));
  genes = colnames(norm.counts)
  inModule = is.finite(match(bwnet$colors,modules))
  modGenes = genes[inModule]
  modTOM=TOM[inModule,inModule]
  dimnames(modTOM)=list(modGenes,modGenes)
  cyt = exportNetworkToCytoscape(modTOM,
                                 edgeFile = paste("merge_CytoscapeInput-edges--1", paste(modules, collapse="-"), ".txt", sep=""),
                                 nodeFile = paste("merge_CytoscapeInput-nodes--1", paste(modules, collapse="-"), ".txt", sep=""),
                                 weighted = TRUE, threshold = -1, nodeNames = modGenes, nodeAttr = bwnet$colors[inModule]);
}


edge_turquoise <- read.delim("merge_CytoscapeInput-edges--1turquoise.txt")
colnames(edge_turquoise)
colnames(edge_turquoise) <- c("source", "target","weight","direction","fromAltName","toAltName")
write.csv(edge_turquoise, file = "edge_turquoise_-1.csv", row.names = FALSE)


node_turquoise <- read.delim("merge_CytoscapeInput-nodes--1turquoise.txt")
colnames(node_turquoise)  
colnames(node_turquoise) <- c("id","altName","node_attributes") 
write.csv(node_turquoise, file = "node_turquoise_-1.csv")


#--------------------------------------------------------------------------------------------------
# Turquoise module filtered network (hub genes are only selected here)
#--------------------------------------------------------------------------------------------------

table(bwnet$colors)

# Filtering edge_turquoise file to retain only rows where both the source and target genes are present in common_hub_genes
filtered_turquoise_edges <- edge_turquoise[
  edge_turquoise$source %in% common_hubgenes_MO & edge_turquoise$target %in% common_hubgenes_MO, ]


filtered_turquoise_edges_0.45 <- filtered_turquoise_edges[filtered_turquoise_edges$weight >= 0.45, ] # Filter edges with weight >= 0.45
write.csv(edge_turquoise, file = "Adult_edge_turquoise_-1(filtered 0.8hubgenes + 0.45weightcut-off).csv", row.names = FALSE)


#-----------------------------------------------------------------------------------
# Functional enrichment analysis (Gene Ontology database)
#------------------------------------------------------------------------------------
 
module.membership.measure <- t(module.membership.measure)
module.membership.measure.pvals <- t(module.membership.measure.pvals)

colnames(module.membership.measure) <- paste("MM", colnames(module_eigengenes), sep = "")
colnames(module.membership.measure.pvals) <- paste("MM.p", colnames(module_eigengenes), sep = "")

colnames(gene.signf.corr.MO) <- paste("GS","V1", sep = "")
colnames(gene.signf.corr.pvals.MO) <- paste("GS.P", "V1", sep = "")

colnames(gene.signf.corr.OPC) <- paste("GS","V1", sep = "")
colnames(gene.signf.corr.pvals.OPC) <- paste("GS.P", "V1", sep = "")


# *****combining the data
allResults_MO = cbind(module.membership.measure, module.membership.measure.pvals, gene.signf.corr.MO, gene.signf.corr.pvals.MO)
head(allResults_MO)

# ******Extract genes in the "turquoise" module
turquoiseModule = "turquoise"
turquoiseGenes = which(bwnet$colors == turquoiseModule)
turquoiseGenes <- as.data.frame(turquoiseGenes)

turquoiseGenes <- unlist(turquoiseGenes)

#is_THAP9_in_turquoise = "THAP9" %in% names(turquoiseGenes)
turquoiseGeneInfo = allResults_MO[turquoiseGenes, ]

turquoise_enrich_BP <- enrichGO(gene = turquoiseGenes,
                                OrgDb = org.Hs.eg.db,
                                ont = "BP",  #BP- biological processes
                                qvalueCutoff = 0.1,
                                readable = TRUE)
### *****View the top enriched terms
head(turquoise_enrich_BP)
dotplot(turquoise_enrich_BP, showCategory = 20) + ggtitle("GO Enrichment Analysis - Biological Process")


turquoise_enrich_MF <- enrichGO(gene = turquoiseGenes,
                                OrgDb = org.Hs.eg.db,
                                ont = "MF",  # MF- molecular function
                                pvalueCutoff = 0.05,
                                qvalueCutoff = 0.1,
                                readable = TRUE)
###   ****View the top enriched terms
head(turquoise_enrich_MF)
dotplot(turquoise_enrich_MF, showCategory = 20) + ggtitle("GO Enrichment Analysis - Molecular Function") 


turquoise_enrich_CC <- enrichGO(gene = turquoiseGenes,
                                OrgDb = org.Hs.eg.db,
                                ont = "CC",  # CC-Cellular component
                                pvalueCutoff = 0.05,
                                qvalueCutoff = 0.1,
                                readable = TRUE)

head(turquoise_enrich_CC)
dotplot(turquoise_enrich_CC, showCategory = 20) + ggtitle("GO Enrichment Analysis - Cellular component")


#----------------------------------------------------------------------------------------------------
# network of selected genes form hub filtered turquoise edge file 
#---------------------------------------------------------------------------------------------------------

# filtered_turquoise_edges is the edge file from turquoise module that containd only the selected hub genes 
genes_of_interest <- c("MBP", "MOG", "THAP9", "SOX6", "SOX4", "SOX5", "SOX11", "PDGFRA")

# genes of interest(goi) *************
Network_goi <- filtered_turquoise_edges[
  filtered_turquoise_edges$source %in% genes_of_interest & filtered_turquoise_edges$target %in% genes_of_interest, ]

write.csv(Network_goi, file ="Network_genes of interest_Adult.csv", row.names = FALSE)


#  *** constructing network from Network_goi file ***  #

# Creating the graph object **************
g <- graph_from_data_frame(Network_goi, directed = FALSE)

edge_label <- E(g)$weight  # Set edge_label to the weights 

node_colors <- brewer.pal(n = vcount(g), name = "Pastel1")   

# Plotting the graph *************
plot(g, 
     vertex.size = 43,            
     vertex.color = node_colors, 
     vertex.label.cex = 0.8,      
     vertex.label.color = "black", 
     vertex.label.font = 2,       
     vertex.shape = "circle",
     edge.width = 2,           
     edge.color = "grey1",     
     #edge.label = round(E(g)$weight, 2), 
     #edge.label.cex = 0.5,        
     #edge.label.color = "black", 
     #edge.label.family = "sans",  
     #edge.label.pos = 0.13,        
     edge.arrow.size = 0.3,          
     edge.arrow.width = 1,
     layout = layout_with_fr(g),  # Fruchterman-Reingold layout
     main = "Network in Adult Samples")  


# saving ***********
save.image(file = "GSE239662_coding+noncoding_Adult.RData")

##############   ******************************************************************************     ########################



