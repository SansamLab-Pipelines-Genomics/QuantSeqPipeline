#!/bin/bash

ml star/2.6.1c

STAR --runThreadN 6 --runMode genomeGenerate --genomeDir starTestIndex --genomeFastaFiles testFiles/test.fa --sjdbGTFfile testFiles/test.gtf --sjdbOverhang 99
