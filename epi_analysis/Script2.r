#####       *********Load libraries ********     ####
library(GenomicRanges)
library(rtracklayer)
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(clusterProfiler)
library(DESeq2)

## -------------------------------------------------------------------------------- ##
# to create Consensus Peaks
## -------------------------------------------------------------------------------- ##

#  ****** Loading replicate peak BED files ********  #
peak_files <- list.files(path = "macs2_output/", pattern = "*.bed", full.names = TRUE)

#  *********Merging peaks across replicates*********  #
peaks <- lapply(peak_files, import)
consensus_peaks <- Reduce(function(x, y) reduce(union(x, y)), peaks)

# Save consensus peaks to BED
export(consensus_peaks, "consensus_peaks.bed")

## -------------------------------------------------------------------------------- ##
# Annotating Peaks 
## -------------------------------------------------------------------------------- ##

peak_file <- "consensus_peaks.bed"

# Loading TxDb for the hg38 genome
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

# To Annotate peaks
peak_anno <- annotatePeak(peak_file, TxDb=txdb, tssRegion=c(-3000, 3000), annoDb="org.Hs.eg.db")
print(peak_anno)
plotAnnoPie(peak_anno)  # Pie chart of genomic feature distribution

###  ***Visualization****  ###

plotAnnoBar(peak_anno) # Genomic distribution of peaks
plotDistToTSS(peak_anno, title="Distribution of Peaks Relative to TSS")

## Saving file *******
write.table(as.data.frame(peak_anno), file="annotated_peaks.txt", sep="\t", row.names=FALSE, quote=FALSE)


## -------------------------------------------------------------------------------- ##
#Differential Acetylation Analysis using DESeq2
## -------------------------------------------------------------------------------- ##

# ***load counts data****  
counts <- read.table("H3K27ac_counts.txt", header = TRUE, row.names = 1)

# ***DESeq2 Analysis***
coldata <- data.frame(condition = factor(c("Adult_OPC", "Adult_MO", "Infant_OPC", "Infant_MO")))
dds <- DESeqDataSetFromMatrix(countData = counts, colData = coldata, design = ~condition)
dds <- DESeq(dds)
res <- results(dds, contrast = c("condition", "Adult_OPC", "Infant_OPC"))

# *** to Save differential peaks results *** #
write.table(res, file = "OligData_differential_peaks.txt", sep = "\t", quote = FALSE)

###############   **************************************************      ##############

