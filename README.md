# Alignment and counting of reads from QuantseqFWD Experiments


- [Scripts](#scripts)
- [Data files in repository](#data-files-in-repository)
- [Additional useful files not included](#additional-useful-files-not-included)
- [Example 1 - run individual scripts on included test files](#example-1---run-individual-scripts-on-included-test-files)
	* [Create a directory for test output](#create-a-directory-for-test-output)
	* [Make test STAR index](#make-test-star-index)
	* [Trim reads](#trim-reads)
	* [Align reads to genome with STAR](#align-reads-to-genome-with-star)
	* [Remove probable internal priming reads and counts over transcripts](#remove-probable-internal-priming-reads-and-counts-over-transcripts)
- [Example 2 - run the entire pipeline on a single sample with slurm wrapper script](#example-2---run-the-entire-pipeline-on-a-single-sample-with-slurm-wrapper-script)

<small><i><a href='http://ecotrust-canada.github.io/markdown-toc/'>Table of contents generated with markdown-toc</a></i></small>


## Scripts
Script | Description
------ | -----------
QuantSeqSlurmWrapper.sh |  Runs all necessary scripts to align and process quantseq reads with slurm.  You must provide a pre-indexed genome for STAR.
TrimReads.sh | Trims reads with bbmap
QuantSeqFilterAndCounts.sh | This will filter reads from a QuantseqFWD experiment that end at a probable genomic poly-A priming site. It will then count the remaining reads over transcript feature coordinates defined in the gtf file.
AlignWithSTAR.sh | Aligns reads to genome with STAR
makeTestSTARIndex.sh | Makes the STAR index for the test data set.
GetInternalPrimingSites.R | This takes a fasta of sequences +/-10bp around putative polyadenylation sites. It checks the genomic sequence in the region [–10..10] surrounding the polyadenylation events and filtered out alignments containing stretches of six consecutive A or with 70% A coverage in any 10-nt sub-window in this region. A bed file of false polyadenylation sites is returned. This is used in QuantSeqFilterAndCounts.sh.

## Data files in repository

File | Description
---- | -----------
**refFiles directory** | 
polyA.fa.gz | fasta with polyA for bbmap trimming
truseq.fa.gz | truseq adapter sequences for bbmap trimming
**testFiles directory** |
test.chrm.lengths | chromosome length for test data set
test.fa | genome sequence for test data set
test.fastq | Quantseq sequencing reads for test data set
test.gtf | transcript feature coordinates for test data set


## Additional useful files not included

####note:  these are not needed for the test run
1. Zebrafish gtf file - from https://www.umassmed.edu/globalassets/lawson-lab/downloadfiles/v4.3.2.gtf
2. Index of Zebrafish genome for STAR - from https://galaxyweb.umassmed.edu/pub/dnext_data/genome_data/zebrafish/GRCz11/refseq_v4.3.2/STARIndex/ 
3. Genome sequence in fasta format

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
### Arguments:
Position | Description | Example for test |
-------- | ----------- | ---------------- |
1 | output directory name | test2_outputdir
2 | name of fastq file with Quantseq reads | testFiles/test.fastq
3 | reference files for trimming | refFiles/polyA.fa.gz,refFiles/truseq.fa.gz
4 | name of directory with the genome index for STAR | starTestIndex/
5 | chromosome sizes file | testFiles/test.chrm.lengths
6 | genome sequence fasta file | testFiles/test.fa
7 | gtf file for transcriptome | testFiles/test.gtf

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
