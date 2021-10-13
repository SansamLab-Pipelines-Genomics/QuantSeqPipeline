#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)
library(zoo)

#' Find likely internal priming sites
#' 
#' This takes a fasta of sequences +/-10bp around putative polyadenylation sites. It checks the genomic sequence in the region [–10..10] surrounding the polyadenylation events and filtered out alignments containing stretches of six consecutive A or with 70% A coverage in any 10-nt sub-window in this region. A bed file of false polyadenylation sites is returned.
#' 
#' @param fastaFile of genomic sequences +/-10bp around the putative polyadenylation sites on the + strand
#' @param outBedFilename name of the .bed file to be written

GetInternalPrimingSites <- function(fastaFile=args[1],
                                    outBedFilename=args[2]){
  flanks <- read.table(fastaFile)
  seqs1 <- flanks[seq(2, nrow(flanks), 2),]
  names(seqs1) <- gsub(">","",flanks[seq(1, nrow(flanks), 2),])
  seqs <- seqs1[-which(duplicated(names(seqs1)))]
  withSixAs <- grep("A{6,}",seqs,perl=T,ignore.case=T)
  set1 <- seqs[withSixAs]
  seqs2 <- seqs[-withSixAs]
  nums <- gsub("A|a",1,seqs2)
  nums <- gsub("[a-zA-Z]",0,nums,perl=T,ignore.case = T)
  nums <- gsub("^0*","",nums,perl=T)
  nums <- gsub("0+$","",nums,perl=T)
  nums_split <- strsplit(nums,split="")
  names(nums_split) <- names(nums)
  nums_split <- lapply(nums_split,as.numeric)
  numAs <- sapply(nums_split,sum)
  nums_split2 <- nums_split[numAs>6]
  lengths <- sapply(nums_split2,length)
  TenOrLess <- nums_split2[lengths<11]
  set2 <- TenOrLess[sapply(TenOrLess,sum)>6]
  nums_split2 <- nums_split2[lengths>10]
  slideSums <- lapply(nums_split2,zoo::rollsum,10)
  slideSumsMax <- sapply(slideSums,max)
  set3 <- (nums_split2)[slideSumsMax>6]
  toRemove <- c(names(set1),names(set2),names(set3))
  seqnames_starts <- as.data.frame(do.call("rbind",strsplit(toRemove,split=":|-",perl=T)))
  seqnames_starts$V3 <- gsub("\\(\\+\\)","",seqnames_starts$V3)
  seqnames_starts$V3 <- gsub("\\(\\-\\)","",seqnames_starts$V3)
  seqnames_starts$V3 <- gsub("\\(","",seqnames_starts$V3)
  seqnames_starts <- data.frame(seqnames_starts$V1,seqnames_starts$V2,seqnames_starts$V3)
  write.table(seqnames_starts,file=outBedFilename,quote=F,sep="\t",col.names=F,row.names=F)
}

GetInternalPrimingSites()