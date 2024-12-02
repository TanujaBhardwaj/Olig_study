###           ********LOAD COUNTS DATA & SAMPLE INFO********        ####

samples <- read.csv("sample_info_Adults_table.csv")
view(samples)

# Read the file
counts.data <- read.csv("Counts_Adults.csv", header = TRUE)
view(counts.data)

counts <- counts.data %>%
  column_to_rownames(var = "GeneID")

dim(counts) # 56852    20
colnames(counts)

counts <- counts[which(rowSums(counts) > 0),] # to remove all the zero counts
dim(counts) #  48949    20

cell_type <- factor(c(samples$cell_type))
coldata <- data.frame(row.names = colnames(counts), cell_type)


####  ************* MAKING dds MATRIX *****************    ###
dds <- DESeqDataSetFromMatrix(countData = counts, colData = coldata, design = ~cell_type)
# Run DESeq2 analysis
dds <- DESeq(dds)


###     ***** Variance stabilizing transformation ******     ####
vsdata <- vst(dds, blind = FALSE)

# ***Plotting PCA*** 
pcaData <- plotPCA(vsdata, intgroup = "cell_type")
plotPCA(vsdata, intgroup = "cell_type")


## ***Dispersion estimates ***
plotDispEsts(dds)

###########   ************* EXTRACTING RESULTS***************   ##################
res <- results(dds, contrast = c("cell_type", "MO", "OPC"))
res
summary(res)
#out of 48948 with nonzero total read count
#adjusted p-value < 0.1
#LFC > 0 (up)       : 8951, 18%
#LFC < 0 (down)     : 10479, 21%
#outliers [1]       : 0, 0%
#low counts [2]     : 15185, 31%
#(mean count < 1)

# log2 fold change (MLE): levels MO vs OPC 
sigs <- na.omit(res)
sigs <- sigs[sigs$padj < 0.05,]
sigs
summary(sigs)
#out of 17711 with nonzero total read count
#adjusted p-value < 0.1
#LFC > 0 (up)       : 8115, 46%
#LFC < 0 (down)     : 9596, 54%
#outliers [1]       : 0, 0%
#low counts [2]     : 0, 0%
#(mean count < 1)

# Save the DESeq2 results 
write.table(sigs, file = "diffExp_MOvsOPC_Adult.txt", sep = ",", row.names = TRUE, quote = FALSE)


# Subset genes with log2FoldChange greater than or equal to 1
sigs_LF1 <- sigs[sigs$log2FoldChange >= 1, ]
dim(sigs_LF1)

###   **************   **VISUALIZTION**   ****************   ##########

plotMA(res) # MA plot

res_deseq <- as.data.frame(res) # Coerce to a data frame
head(res_deseq)

res_deseq$significant <- ifelse(res_deseq$padj < .1, "Significant", NA)

## *** detailed MA Plot ****
ggplot(res_deseq, aes(baseMean, log2FoldChange, colour=padj)) + 
  geom_point(size=1) + scale_y_continuous(limits=c(-3, 3), oob=squish) + 
  scale_x_log10() + geom_hline(yintercept = 0, colour="darkorchid4", size=1, linetype="longdash") + 
  labs(x="mean of normalized counts", y="log fold change",                            # Y-axis label
       title="MA plot Adult samples") + scale_colour_viridis(direction=-1, trans='sqrt') + 
  theme_bw() + 
  geom_density_2d(colour="black", size=1) +
  # Add a cross marker and label for "THAP9"
  geom_point(data = subset(res_deseq, rownames(res_deseq) == "THAP9"),
             aes(baseMean, log2FoldChange), shape = 19, size = 2, color = "red") +  
  geom_text(data = subset(res_deseq, rownames(res_deseq) == "THAP9"),
            aes(baseMean, log2FoldChange, label = "THAP9"), vjust = -1, size = 5, color = "red") 

######        ****FOR Ploting counts for GENE****       ######

THAPCounts <- plotCounts(dds, gene="THAP9", intgroup="cell_type", returnData=TRUE)  # Extract counts for the gene THAP9

# Plot the data using ggplot2
colourPallette <- c("#7145cd","#90de4a","#cd46c1","#77dd8e","#592b79","#d7c847","#6378c9","#619a3c",
                    "#d44473","#63cfb6","#dd5d36","#5db2ce","#8d3b28","#b1a4cb","#af8439","#c679c0","#4e703f",
                    "#753148","#cac88e","#352b48","#cd8d88","#463d25","#556f73","#bbcfc4")

ggplot(THAPCounts, aes(x=cell_type, y=count, color=cell_type, group=cell_type)) + 
  geom_point(size=3) +                               
  geom_line() +                                    
  theme_bw() + 
  theme(axis.text.x=element_text(angle=15, hjust=1)) + 
  scale_colour_manual(values=colourPallette) +       
  guides(colour=guide_legend(ncol=3)) +              
  ggtitle("THAP9 Counts by Cell Type Adult samples") +          
  ylab("Counts") +                                   
  xlab("Cell Type")                                 


###########     ******Plotting results for THAP9 gene*****      ##########

thap9_res <- subset(res, rownames(res) == "THAP9")  # Filter for THAP9 gene in the DESeq2 results
sampleData <- samples

colData(vsdata)$cell_type <- sampleData$cell_type[match(colnames(vsdata), sampleData$sample_id)]
colnames(colData(vsdata))

thap9_expr <- assay(vsdata)["THAP9", ]  
df_thap9 <- data.frame(Expression = thap9_expr,
                       Condition = colData(vsdata)$cell_type)

head(df_thap9)

###      ***Create the box plot***
ggplot(df_thap9, aes(x = cell_type, y = Expression, fill = cell_type)) +
  geom_boxplot() +
  labs(title = "Expression of THAP9 in MOs and OPCs Adult samples",
       x = "Sample", y = "Expression") +
  theme_minimal() +
  scale_fill_brewer(palette = "Set2") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


################      *****PCA plots******           ##################
percentVar <- round(100 * attr(pcaData, "percentVar"))
pcaData_df <- pcaData$data
str(pcaData_df)

pdf('pca_plot_Adult.pdf')
ggplot(pcaData_df, aes(x = PC1, y = PC2, color = cell_type, label = name)) +
  geom_point(size = 3) +
  geom_text(vjust = 1.5, size = 3.5) + # This adds text labels to the points
  xlab(paste0("PC1 ")) +
  ylab(paste0("PC2 ")) +
  labs(title = "PCA Plot of Adult samples") +
  theme_bw()
dev.off()

############################################################# 

###   ***Create the volcano plot with annotation ***  ###
volcano_Adult <- ggplot(sigs) +
  aes(y = -log10(padj), x = log2FoldChange) +
  geom_point(size = 2) +
  geom_hline(yintercept = -log10(0.01), linetype = "longdash", colour = "grey", size = 1) +
  geom_vline(xintercept = 1, linetype = "longdash", colour = "#BE6", size = 1) +
  geom_vline(xintercept = -1, linetype = "longdash", colour = "#BE684D", size = 1) +
  annotate("rect", xmin = 1, xmax = Inf, ymin = -log10(0.01), ymax = Inf, alpha = .2, fill = "#BE6") +
  annotate("rect", xmin = -Inf, xmax = -1, ymin = -log10(0.01), ymax = Inf, alpha = .2, fill = "#BE684D") +
  labs(title = "Volcano plot",
       subtitle = "MOs v/s OPCs in Adult samples") +
  theme_bw() +
  geom_point(data = subset(sigs, rownames(sigs) == "THAP9"), 
             shape = 8, size = 2, color = "tomato") + 
  geom_text(data = subset(sigs, rownames(sigs) == "THAP9"),
            aes(label = "THAP9"), vjust = -1, color = "tomato")


print(volcano_Adult)

##########    ******* To make the expression graph ***********   ###########
df_thap9 <- data.frame(
  Expression = thap9_expr,
  Condition = cell_type
)

# Extractting  statistical results
log2FC <- round(thap9_res$log2FoldChange, 2)
p_adj <- format(thap9_res$padj, scientific = TRUE)

# ***Plotting ***
ggplot(df_thap9, aes(x = Condition, y = Expression, fill = Condition)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.9, width = 0.5, color = "black", size = 0.8) +
  labs(
    title = "Expression of THAP9 in MOs vs OPCs (Adult)",
    subtitle = paste("Log2 Fold Change = ", log2FC, 
                     "\nAdjusted p-value = ", p_adj),
    x = "Sample",
    y = "Expression Level"
  ) +
  scale_fill_manual(
    values = c("MO" = "#66CCFF", "OPC" = "#FF6666"),
    labels = c("OPCs", "MOs")
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    plot.subtitle = element_text(hjust = 0.5, size = 12, margin = margin(b = 10)),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 12),
    legend.title = element_blank(),
    legend.position = "top",  
    legend.text = element_text(size = 12),
    panel.grid.minor = element_blank(),  
    panel.grid.major.x = element_blank(),  
    panel.border = element_rect(fill = NA, color = "gray80", size = 0.5),  
    plot.background = element_rect(fill = "#F5F5F5", color = NA)  
  ) +
  coord_cartesian(ylim = c(7, 10))

############   ***********************************     ###############

# upregulated and downregulated genes
up_genes <- rownames(subset(res, log2FoldChange > 0 & padj < 0.1))
down_genes <- rownames(subset(res, log2FoldChange < 0 & padj < 0.1))

# THAP9 status within upregulated or downregulated genes
thap9_status <- ifelse("THAP9" %in% up_genes, "Upregulated", 
                       ifelse("THAP9" %in% down_genes, "Downregulated", "Not found"))


fill_colors <- c("green4", "salmon3")  # Upregulated in green, Downregulated in salmon
highlight_color <- "red"  # Highlight color for THAP9

# ******** Venn diagram with two sets: upregulated and downregulated *********** #
venn.plot <- venn.diagram(
  x = list(Upregulated = up_genes, Downregulated = down_genes),
  category.names = c("Upregulated", "Downregulated"),
  filename = NULL,
  fill = fill_colors,
  alpha = 0.6,
  cex = 1.8,
  fontface = "bold",
  cat.cex = 1.5,
  cat.fontface = "bold",
  cat.pos = c(-5, 5),  
  cat.dist = c(0.05, 0.05)  
)
grid.newpage()
grid.draw(venn.plot)

if (thap9_status == "Upregulated") {
  grid.text("THAP9", x = unit(0.76, "npc"), y = unit(0.55, "npc"), 
            gp = gpar(col = highlight_color, fontsize = 13, fontface = "bold"))
} else if (thap9_status == "Downregulated") {
  grid.text("THAP9", x = unit(0.7, "npc"), y = unit(0.7, "npc"), 
            gp = gpar(col = highlight_color, fontsize = 13, fontface = "bold"))
} else {
  message("THAP9 not found in any category.")
}

grid.text("Adult Samples (MOs vs OPCs)", 
          x = unit(0.5, "npc"), y = unit(0.95, "npc"), 
          gp = gpar(col = "black", fontsize = 16, fontface = "bold"))
grid.text("log2FoldChange > 0 & padj < 0.1", 
          x = unit(0.5, "npc"), y = unit(0.05, "npc"), 
          gp = gpar(col = "black", fontsize = 12, fontface = "italic"))


#####################               ***************************************              #########################
