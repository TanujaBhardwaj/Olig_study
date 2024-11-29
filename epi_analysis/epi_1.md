##### ********************* SAMTOOLS  ************************************ ####

# to convert the .sam files to .bam files using samtools
	samtools view -S -b ##.sam > ##.bam


#### ************** PEAKCALLING MACS2 ************************************ ####

# Using MACS2 for peak calling

	macs2 callpeak -t chip.bam \
    -c control.bam \
    -f BAM -g hs \
    -n chip_Sample \
    --outdir macs2 2> macs2.log


# ****using bedtools to remove the blacklist regions*****  #
	bedtools intersect \
    -v \
    -a macs2/sample_peaks.narrowPeak \
    -b ../reference_data/blacklist.v2.bed \
    macs2/psample_peaks_filtered.bed

