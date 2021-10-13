#!/bin/bash

ml star/2.6.1c

# Aligns reads to genome with STAR
#
# $1 - name of fastq file with Quantseq reads
# $2 - directory with indexed genome
# $3 - output directory name (must exist already)
#
# Example:
#
# AlignWithStar.sh \
# test.fastq \
# starTestIndex/ \
# test_outputdir

FASTQ=$1
GENOME=$2
OUT=$(basename ${FASTQ} .fastq)
OUTPUTDIR=$3

STAR \
--genomeDir ${GENOME} \
--runThreadN 6 \
--readFilesIn ${FASTQ} \
--outFileNamePrefix ${OUTPUTDIR}/${OUT} \
--outSAMtype BAM SortedByCoordinate \
--outSAMunmapped Within \
--outSAMattributes Standard \
--alignIntronMin 20 \
--alignIntronMax 1000000 \
