##############    ************** Downloading and indexing***********      #####################

# Downloading files using sratools
	
	fastq-dump  
	vdb-config --prefetch-to-cwd
	prefetch SRR#####
	fasterq-dump SR#### --split-files --skip-technical


######  -------------------------------------------------------------------------    #########

#### Using fastp for trimming and removing low quality reads **************

for paired end files:
	fastp -i SRR####_1.fastq.gz -I SRR####_2.fastq.qz -o ####_1.fastq.gz -O ####_2.fastq.qz -l 36 --detect_adapter_for_pe -c


##############################################################################
### Quality control *************
	multiqc .fastq*


######  -------------------------------------------------------------------------    #########

###############**************************Indexing**************************************###################

# Using Kallisto ***********:
	kallisto index -i Homo_sapiens.GRCh38.indexed.fa Homo_sapiens.GRCh38.dna.alt 


## Using Bowtie ****************
	bowtie2-build hg38.fa indexed_hg38


# bowtie2 alignment:
	bowtie2 -x reference_genome -1 SRR#####_1.fastq -2 SRR####_2.fastq -S SRR#####.sam