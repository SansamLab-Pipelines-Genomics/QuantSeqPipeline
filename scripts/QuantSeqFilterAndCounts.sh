#!/bin/bash

ml samtools
ml R 
ml bedtools
ml bedops
ml picard
ml python
ml htseq

# Filters aligned reads and counts over transcripts.
#
# This will filter reads from a QuantseqFWD experiment that end at a probable genomic poly-A priming site. It will then count the remaining reads over gtf transcripts. Takes four arguments.
#
# $1 - basename of .bam file with aligned Quantseq reads
# $2 - chromosome sizes file
# $3 - genome sequence fasta file
# $4 - gtf file for transcriptome
# $5 - directory where outputs should be saved
#
# Examples:
#
# QuantSeqFilterAndCounts.sh \
# test \
# testFiles/test.chrm.lengths \
# testFiles/test.fa \
# testFiles/test.gtf \
# test_outputdir
#
# QuantSeqFilterAndCounts.sh \
# 1_WT_1-cell_S1 \
# /Volumes/shared-refs/star_genomes/grcz11-zebrafish/refseq_v4.3.2/chrNameLength.txt \
# /Volumes/shared-refs/zebrafish/danRer11/danRer11.fa \
# v4.3.2.gtf \
# 1_WT_1-cell_S1_outputdir

BSENME=$1
CHROMSIZES=$2
GENOMEFASTA=$3
GTFFILE=$4
OUTPUTDIR=$5

# check if script is started via SLURM or bash
# if with SLURM: there variable '$SLURM_JOB_ID' will exist
# `if [ -n $SLURM_JOB_ID ]` checks if $SLURM_JOB_ID is not an empty string
if [ -n $SLURM_JOB_ID ];  then
    # check the original location through scontrol and $SLURM_JOB_ID
    SCRIPT_PATH=$(scontrol show job $SLURM_JOBID | awk -F= '/Command=/{print $2}')
else
    # otherwise: started with bash. Get the real location.
    SCRIPT_PATH=$(realpath $0)
fi

# getting location of script
SCRIPTDIR=${SCRIPT_PATH%%[[:space:]]*}
SCRIPTDIR=$(dirname "${SCRIPTDIR}")

echo "################"
echo $SCRIPT_PATH
echo ${SCRIPTDIR}
echo "################"

#SCRIPTDIR=$(dirname "${BASH_SOURCE[0]}")

main () {
	SubsetBamByStrand
	ConvertStrandBamsToBeds
	GetSeqsAtReadsEnds
	GetReadsPrimedAtAs
	MakeBamsWithoutReadsPrimedAtAs
	CountReadsOverTranscripts
}


SubsetBamByStrand () {
	# get positive strand reads bam BSENME
	samtools view -F 16 -b \
	-o ${OUTPUTDIR}/${BSENME}_positive_strand.bam \
	${OUTPUTDIR}/${BSENME}Aligned.sortedByCoord.out.bam

	# get negative strand reads bam BSENME
	samtools view -f 16 -b \
	-o ${OUTPUTDIR}/${BSENME}_negative_strand.bam \
	${OUTPUTDIR}/${BSENME}Aligned.sortedByCoord.out.bam

	# index positive strand reads bam BSENME
	samtools index ${OUTPUTDIR}/${BSENME}_positive_strand.bam

	# index negative strand reads bam BSENME
	samtools index ${OUTPUTDIR}/${BSENME}_negative_strand.bam
}


ConvertStrandBamsToBeds () {
	# convert positive strand reads bam to bed +/- 10bp around ends
	bedtools bamtobed -i ${OUTPUTDIR}/${BSENME}_positive_strand.bam |\
	bedtools flank -i stdin \
	-g ${CHROMSIZES} \
	-l 0 -r 10 |
	bedops --range -10:0 --everything - > ${OUTPUTDIR}/${BSENME}_positive_strand_flank.bed

	# convert negative strand reads bam to bed +/- 10bp around starts
	bedtools bamtobed -i ${OUTPUTDIR}/${BSENME}_negative_strand.bam |\
	bedtools flank -i stdin \
	-g ${CHROMSIZES} \
	-l 0 -r 10 |
	bedops --range -10:0 --everything - > ${OUTPUTDIR}/${BSENME}_negative_strand_flank.bed
}
	

GetSeqsAtReadsEnds () {
	# make fasta of genomic sequence +/-10bp around end 
	bedtools getfasta -s -fi ${GENOMEFASTA} \
	-bed ${OUTPUTDIR}/${BSENME}_positive_strand_flank.bed > ${OUTPUTDIR}/${BSENME}_positive_strand_flank.fasta

	# make fasta of genomic sequence +/-10bp around end 
	bedtools getfasta -s -fi ${GENOMEFASTA} \
	-bed ${OUTPUTDIR}/${BSENME}_negative_strand_flank.bed > ${OUTPUTDIR}/${BSENME}_negative_strand_flank.fasta
}


GetReadsPrimedAtAs () {
	# id likely internal priming sites for positive strand reads
	Rscript ${SCRIPTDIR}/GetInternalPrimingSites.R \
	${OUTPUTDIR}/${BSENME}_positive_strand_flank.fasta \
	${OUTPUTDIR}/${BSENME}_positive_strand_internalPriming.bed

	# id likely internal priming sites for negative strand reads
	Rscript ${SCRIPTDIR}/GetInternalPrimingSites.R \
	${OUTPUTDIR}/${BSENME}_negative_strand_flank.fasta \
	${OUTPUTDIR}/${BSENME}_negative_strand_internalPriming.bed

	# get read names at internal priming sites
	bedtools intersect -wa -f 1 \
	-a ${OUTPUTDIR}/${BSENME}_positive_strand_flank.bed \
	-b ${OUTPUTDIR}/${BSENME}_positive_strand_internalPriming.bed | \
	cut -f4 > ${OUTPUTDIR}/${BSENME}_positive_strand_internalPrimingReads.txt
	bedtools intersect -wa -f 1 \
	-a ${OUTPUTDIR}/${BSENME}_negative_strand_flank.bed \
	-b ${OUTPUTDIR}/${BSENME}_negative_strand_internalPriming.bed | \
	cut -f4 > ${OUTPUTDIR}/${BSENME}_negative_strand_internalPrimingReads.txt
}


MakeBamsWithoutReadsPrimedAtAs () {
	# extract reads NOT at internal priming sites - positive strand
	picard FilterSamReads \
	I=${OUTPUTDIR}/${BSENME}_positive_strand.bam \
	O=${OUTPUTDIR}/${BSENME}_positive_strand_filtered.bam \
	READ_LIST_FILE=${OUTPUTDIR}/${BSENME}_positive_strand_internalPrimingReads.txt \
	FILTER=excludeReadList
	# extract reads NOT at internal priming sites - negative strand
	picard FilterSamReads \
	I=${OUTPUTDIR}/${BSENME}_negative_strand.bam \
	O=${OUTPUTDIR}/${BSENME}_negative_strand_filtered.bam \
	READ_LIST_FILE=${OUTPUTDIR}/${BSENME}_negative_strand_internalPrimingReads.txt \
	FILTER=excludeReadList

	# recombine positive and negative bams
	samtools merge -f --write-index \
	${OUTPUTDIR}/${BSENME}_filtered.bam \
	${OUTPUTDIR}/${BSENME}_positive_strand_filtered.bam \
	${OUTPUTDIR}/${BSENME}_negative_strand_filtered.bam
}

CountReadsOverTranscripts () {
	htseq-count -f bam -s yes ${OUTPUTDIR}/${BSENME}_filtered.bam ${GTFFILE} > ${OUTPUTDIR}/${BSENME}_counts.tab
}

main