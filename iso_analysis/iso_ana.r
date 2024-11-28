##############     *** Load libraries ***    ########
library(tidyverse)
library(datapasta)
library(IsoformSwitchAnalyzeR)
library(tximport)
library(rhdf5)
library(dplyr)
library(cowplot)

###    ****load sample files****   ###
targets <- read.csv("Infant&Adult_sample.txt", sep = "\t")

targets.mod <- targets %>%
  dplyr::select(1,6) %>%
  dplyr::rename(sampleID = sample, condition = cell_type) %>%
  dplyr::select(sampleID, condition)

dim(targets.mod)

##   **** import data  ****  ##
path <- file.path(targets.mod$sampleID, "abundance.tsv") 
all(file.exists(path))
which(file.exists(path))

Txi_trans <- importIsoformExpression(sampleVector = path)
dim(Txi_trans$counts)

sampleLabels <- targets.mod$sampleID

colnames(Txi_trans$abundance) <- c("isoform_id", sampleLabels)
colnames(Txi_trans$counts) <- c("isoform_id", sampleLabels)


# *****to make a switchlist for analysis****
mySwitchList <- importRdata(
  isoformCountMatrix   = Txi_trans$counts,
  isoformRepExpression = Txi_trans$abundance,
  designMatrix         = targets.mod,
  removeNonConvensionalChr = TRUE,
  addAnnotatedORFs=TRUE,
  ignoreAfterPeriod=TRUE,
  isoformExonAnnoation = "Homo_sapiens.GRCh38.112.chr_patch_hapl_scaff.gtf.gz",
  isoformNtFasta = "Homo_sapiens.GRCh38.cdna.all.fa.gz",
  fixStringTieAnnotationProblem = TRUE,
  showProgress = TRUE)

summary(mySwitchList)

head(mySwitchList$conditions,2)
# condition nrReplicates
#1   MO_Adult            3
#2  MO_Infant            4
#3  OPC_Adult            3
#4 OPC_Infant            4

names(mySwitchList)

saveRDS(mySwitchList, file = "mySwitchList_adult_&_infants.RDS")



#  *****Pre-filter to remove all the lowely expressed features******** 
mySwitchList <- preFilter(
  switchAnalyzeRlist = mySwitchList,
  geneExpressionCutoff = 1,
  isoformExpressionCutoff = 0,
  removeSingleIsoformGenes = TRUE
)

####   *********************** DEXSeq ***********************    ####
mySwitchList_DEXSeq_analysed <- isoformSwitchTestDEXSeq(
  switchAnalyzeRlist = mySwitchList,
  reduceToSwitchingGenes=FALSE)

extractSwitchSummary(mySwitchList_DEXSeq_analysed)

######################   *****************     #################

###    *****Switchplot in Infant samples*******    ####
switchPlot(
  mySwitchList_DEXSeq_analysed,
  gene = 'THAP9',
  condition1 = 'MO_Infant',
  condition2 = "OPC_Infant",
  localTheme = theme_bw()
)


###     *****Switchplot in Adult samples*******   ######
switchPlot(
  mySwitchList_DEXSeq_analysed,
  gene = 'THAP9',
  condition1 = 'MO_Adult',
  condition2 = "OPC_Adult",
  localTheme = theme_bw()
)

####   ******** individual plots ********    ####
switchPlotTranscript(mySwitchList_DEXSeq_analysed, gene = 'THAP9')
switchPlotGeneExp(mySwitchList_DEXSeq_analysed, gene = 'THAP9',
                  condition1 = 'MO_Adult',
                  condition2 = "OPC_Adult",
                  localTheme = theme_bw())

switchPlotIsoExp(mySwitchList_DEXSeq_analysed, gene = 'THAP9',
                 condition1 = 'MO_Infant',
                 condition2 = "OPC_Infant",
                 localTheme = theme_bw())

switchPlotIsoUsage(mySwitchList_DEXSeq_analysed, gene = 'THAP9',
                   condition1 = 'MO_Infant',
                   condition2 = "OPC_Infant",
                   localTheme = theme_bw())

###########   *******************external analysis***************************    ######################

#### for external analysis we made a short file of only the genes we are interested in
# to annotate the isoforms involved in the isoform switches and the annotation of Open Reading Frames (ORFs).

# because external analysis tool does not tend to support heavy files we will subset data
genes_of_interest <- c(
  'THAP9', 'THAP12', 'OLIG2', 'SOX10', 'THAP6', 
  'THAP1', 'THAP2', 'THAP3', 'THAP4', 'THAP5', 
  'THAP7', 'THAP8', 'ZNF24', 'POU5F1'
)

### Subsetting  *******
mySwitchList_filt <- subsetSwitchAnalyzeRlist(
  mySwitchList_DEXSeq_analysed,
  mySwitchList_DEXSeq_analysed$isoformFeatures$gene_name %in% genes_of_interest
)

head(mySwitchList_filt$isoformFeatures)

## To check ORFs
"orfAnalysis" %in% names(mySwitchList_filt)


###    ******extract sequences*******    ####
mySwitchList_filt <- extractSequence(
  mySwitchList_filt, 
  pathToOutput = 'sequences/',
  writeToFile=TRUE # to avoid output when running this example data
)

###     **to Add CPC2 analysis**     ####
mySwitchList_filt_analyzed <- analyzeCPC2(
  switchAnalyzeRlist   = mySwitchList_filt,
  pathToCPC2resultFile = "result_cpc2_Oligo.txt",
  removeNoncodinORFs   = TRUE  
)

###     **to Add PFAM analysis**     ####
mySwitchList_filt_analyzed <- analyzePFAM(
  switchAnalyzeRlist   = mySwitchList_filt_analyzed,
  pathToPFAMresultFile = "Pfam_Oligo_results.txt",
  showProgress=FALSE
)

###      **to Add SignalP analysis**   ####
mySwitchList_filt_analyzed <- analyzeSignalP(
  switchAnalyzeRlist       = mySwitchList_filt_analyzed,
  pathToSignalPresultFile  = "prediction_results_SignalP_Oligo.txt"
)

###       **to Add IUPred2A analysis**  ####
mySwitchList_filt_analyzed <- analyzeIUPred2A(
  switchAnalyzeRlist        = mySwitchList_filt_analyzed,
  pathToIUPred2AresultFile = "IUPred_Oligo.result",
  showProgress = FALSE
)

###       **to Add DeepLoc2 analysis**  ####
mySwitchList_filt_analyzed <- analyzeDeepLoc2(
  switchAnalyzeRlist = mySwitchList_filt_analyzed,
  pathToDeepLoc2resultFile = "results_DeepLoc2_Oligo.csv",
  quiet = FALSE
)

###       **to Add DeepTMHMM analysis**  ####
mySwitchList_filt_analyzed <- analyzeDeepTMHMM(
  switchAnalyzeRlist   = mySwitchList_filt_analyzed,
  pathToDeepTMHMMresultFile = "DeepTMHMM_Oligo_results.gff3",
  showProgress=FALSE
)

mySwitchList_filt_analyzed
#This switchAnalyzeRlist list contains:
#82 isoforms from 14 genes
#6 comparison from 4 conditions (in total 14 samples)
#Feature analyzed:
#  [1] "Isoform Switch Identification, ORFs, ntSequence, aaSequence, 
#Protein Domains, IDR, Sub-cellular localization, topologyAnalysis, Coding Potential"

mySwitchList_filt_analyzed <- analyzeAlternativeSplicing(
  switchAnalyzeRlist = mySwitchList_filt_analyzed,
  quiet=TRUE
)

table( mySwitchList_filt_analyzed$AlternativeSplicingAnalysis$IR ) #number of intron retentions (IR)
consequencesOfInterest <- c('intron_retention','coding_potential','NMD_status','domains_identified','ORF_seq_similarity')


###    *** to analyse the switch consequences ***  ###
mySwitchList_filt_analyzed <- analyzeSwitchConsequences(
  mySwitchList_filt_analyzed,
  consequencesToAnalyze = consequencesOfInterest, 
  dIFcutoff = 0.1, 
  showProgress=FALSE
)

extractSwitchSummary(mySwitchList_filt_analyzed, dIFcutoff = 0.1, filterForConsequences = FALSE)

###   **** to extract the splice summary ****  ###
extractSplicingSummary(
  mySwitchList_filt_analyzed,
  asFractionTotal = FALSE,
  plotGenes=FALSE
)

#### ********splicingEnrichment******** #########
splicingEnrichment <- extractSplicingEnrichment(
  mySwitchList_filt_analyzed,
  splicingToAnalyze='all',
  returnResult=TRUE,
  returnSummary=TRUE
)

extractSplicingEnrichmentComparison(
  mySwitchList_filt_analyzed,
  splicingToAnalyze=c('A5','MES','ATSS'), 
  returnResult=FALSE 
)

extractSplicingGenomeWide(
  mySwitchList_filt_analyzed,
  featureToExtract = 'all',  # all isoforms stored in the switchAnalyzeRlist
  splicingToAnalyze = c('A5','MES','ATTS'), 
  plot=TRUE,
  returnResult=FALSE  
)


##########   ********* final plots ********    ##########

# Generate and save the plot *****
png("./Infant_analysed.png", 
    width = 1200, height = 800)

switchPlot(
  mySwitchList_filt_analyzed,
  gene = 'THAP9',
  condition1 = 'MO_Infant',
  condition2 = "OPC_Infant",
  localTheme = theme_bw()
)
dev.off()


png("./Adult_analysed.png", 
    width = 1200, height = 800)

switchPlot(
  mySwitchList_filt_analyzed,
  gene = 'THAP9',
  condition1 = 'MO_Adult',
  condition2 = "OPC_Adult",
  localTheme = theme_bw()
)
dev.off()

###############         ***************************************************    ###########################


