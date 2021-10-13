# Alignment and counting of reads from QuantseqFWD Experiments

## Scripts
1. AlignWithSTAR.sh
2. GetInternalPrimingSites.R
3. makeTestSTARIndex.sh
4. QuantSeqFilterAndCounts.sh

## Data files in repository
1. chromosome.sizes
2. small fastq file
3. expected count results
4. Zebrafish gtf file - from https://www.umassmed.edu/globalassets/lawson-lab/downloadfiles/v4.3.2.gtf; 

## Additional data files required
1. genome sequence in fasta format
2. genome sequence indexed for STAR


## Example 1 - run individual scripts on included test files

### Create a directory for test output
```
mkdir test_outputdir
```

### Make test STAR index
```
mkdir starTestIndex
sbatch --cpus-per-task 6 --mem 48G scripts/makeTestSTARIndex.sh
```

### Trim reads
```
sbatch scripts/TrimReads.sh \
testFiles/test.fastq \
starTestIndex/ \
test_outputdir \
refFiles/polyA.fa.gz,refFiles/truseq.fa.gz
```

### Align reads to genome with STAR
```
sbatch --cpus-per-task 6 --mem 48G scripts/AlignWithSTAR.sh \
testFiles/test.fastq \
starTestIndex/ \
test_outputdir
```

### Remove probable internal priming reads and counts over transcripts
```
sbatch scripts/QuantSeqFilterAndCounts.sh \
test \
testFiles/test.chrm.lengths \
testFiles/test.fa \
testFiles/test.gtf \
test_outputdir
```

## Example 2 - run the entire pipeline on a single sample with slurm wrapper script
### Arguments
1. output directory name (test2_outputdir)
2. name of fastq file with Quantseq reads (testFiles/test.fastq)
3. reference files for trimming (refFiles/polyA.fa.gz,refFiles/truseq.fa.gz)
4. name of directory with the genome index for STAR (starTestIndex/)
5. chromosome sizes file (testFiles/test.chrm.lengths)
6. genome sequence fasta file (testFiles/test.fa)
7. gtf file for transcriptome (testFiles/test.gtf)


```
sbatch scripts/QuantSeqSlurmWrapper.sh \
test2_outputdir \
testFiles/test.fastq \
refFiles/polyA.fa.gz,refFiles/truseq.fa.gz \
starTestIndex/ \
testFiles/test.chrm.lengths \
testFiles/test.fa \
testFiles/test.gtf
```
