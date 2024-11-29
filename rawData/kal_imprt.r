# __________________________________________________________________________________________________  # 
             ##### ********* IMPORTING DATA FOR ANALYSIS *************** #####
# ___________________________________________________________________________________________________ #

#   **********Load packages********* ----
library(rhdf5) 
library(tidyverse) 
library(tximport)
library(ensembldb) 
library(EnsDb.Hsapiens.v86) 


###### Study design  ********** ----
targets <- read_tsv("sample_Info.txt")

path <- file.path(targets$sample, "abundance.tsv") # setting file paths to mapped data

all(file.exist(path)) # to check if the path is correct 

# ******TO Get annotation information********* ----
Tx <- transcripts(EnsDb.Hsapiens.v86, columns=c("tx_id", "gene_name"))
Tx <- as_tibble(Tx)   # to change it into tabular form
Tx <- dplyr::rename(Tx, target_id = tx_id)
Tx <- dplyr::select(Tx, "target_id", "gene_name") 


# ****importing Kallisto counts into R using Tximport***** ----
Txi_gene <- tximport(path,
	type = "Kallisto",
	tx2gene = Tx,
	txOut = FALSE,
	countsFromAbundance = "lengthScaledTPM",
	ignoreTxVersion = TRUE)
beep(sound = 3)

class(Txi_gene)
names(Txi_gene)

sampleLabels <- targets$SampleID
Counts_data <- Txi_gene$counts
colSums(Counts_data)
colnames(Counts_data) <- sampleLabels

## *****************  ###
Counts_data <- as_tibble(myCounts, rownames = "geneID")
Counts_data
colnames(Counts_data) <- c("geneID", sampleLabels)

####   ****** Saving file******    ####
write.csv(Counts_data, "Counts_data.csv")


#################     *************************************************************      #########################
