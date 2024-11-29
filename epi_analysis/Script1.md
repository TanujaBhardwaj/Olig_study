# Script for ChIPseq H3K27ac processing

## Indexing using bowtie2
    bowtie2-build hg38.fa hg38

## Align reads using bowtie2 
    bowtie2 -x hg38 -U file.fastq.gz -S aligned.sam -p 8

## Converting SAM to BAM using samtools
    samtools view -bS aligned.sam | samtools sort -o aligned_sorted.bam
    rm aligned.sam  # to remove SAM file

# Filtering
    samtools view -F 2048 -f 2 -q 30 -b aligned_sorted.bam > filtered.bam

# To remove ENCODE blacklisted regions using bedtools
    bedtools intersect -v -abam filtered.bam -b encode_blacklist.bed > filtered_clean.bam

# MACS2 for peak calling
    macs2 callpeak -B -f BAMPE --broad --broad-cutoff 0.00001 -t filtered_clean.bam -c input_control.bam -n sample --outdir macs2_output

# To generate the read counts for peaks
#### To combine BAM files and calculate coverage
    bedtools coverage -a consensus_peaks.bed -b filtered_clean.bam > read_counts.txt

# To convert BAM to BigWig file using deepTools
    bamCoverage -b filtered_clean.bam -o sample.bw --normalizeUsing RPGC --effectiveGenomeSize 2913022398 --binSize 10
