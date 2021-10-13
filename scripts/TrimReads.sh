#!/bin/bash

ml bbmap

# Aligns reads to genome with STAR
#
# $1 - name of fastq file with Quantseq reads
# $2 - output directory name (must exist already)
# $3 - reference files
#
# Example:
#
# TrimReads.sh \
# test.fastq \
# test_outputdir \
# refFiles/polyA.fa.gz,refFiles/truseq.fa.gz
# 
# TrimReads.sh \
# 2_WT_1-cell_S2_R1_001.fastq \
# starTestIndex/ \
# test_outputdir
# "polyA.fa.gz,truseq.fa.gz" 


FASTQ=$1
REFERENCEFILES=$3
BASENAME=$(basename ${FASTQ} .fastq)
OUTPUTDIR=$2

bbduk.sh in=${FASTQ} out=${OUTPUTDIR}/${BASENAME}_trimmed.fastq ref=${REFERENCEFILES}
