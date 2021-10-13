#!/bin/bash

# Aligns reads to genome with STAR
#
# $1 - output directory name
# $2 - name of fastq file with Quantseq reads
# $3 - reference files for trimming
# $4 - name of directory with the genome index for STAR
# $5 - chromosome sizes file
# $6 - genome sequence fasta file
# $7 - gtf file for transcriptome
#
# Example:
#
# QuantSeqSlurmWrapper.sh \
# test2_outputdir \
# testFiles/test.fastq \
# refFiles/polyA.fa.gz,refFiles/truseq.fa.gz \
# starTestIndex/ \
# testFiles/test.chrm.lengths \
# testFiles/test.fa \
# testFiles/test.gtf

OUTPUTDIR=$1
FASTQ=$2
TRIMMINGREFS=$3
STARINDEXDIR=$4
CHROMSIZES=$5
GENOMEFASTA=$6
GTFFILE=$7

BASENAME=$(basename ${FASTQ} .fastq)

mkdir -p ${OUTPUTDIR}

# Trim reads
jid1=$(\
sbatch \
--parsable \
scripts/TrimReads.sh \
${FASTQ} \
${OUTPUTDIR} \
${TRIMMINGREFS})

# Align reads to genome with STAR
jid2=$(\
sbatch \
--parsable \
--dependency=afterany:${jid1} \
--cpus-per-task 6 \
--mem 48G \
scripts/AlignWithSTAR.sh \
${OUTPUTDIR}/${BASENAME}_trimmed.fastq \
${STARINDEXDIR} \
${OUTPUTDIR})

# Remove probable internal priming reads and counts over transcripts
sbatch \
--dependency=afterany:${jid2} \
scripts/QuantSeqFilterAndCounts.sh \
${BASENAME}_trimmed \
${CHROMSIZES} \
${GENOMEFASTA} \
${GTFFILE} \
${OUTPUTDIR}
