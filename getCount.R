library("Rsamtools")
library(IRanges)
library(GenomicRanges)

setwd('/scratch/cqs/guoy1/cleveland/count')

args<-commandArgs(trailingOnly=TRUE)
bamfile<-args[1]
sample<-args[2]

#Return GRangesList object in both gene level and transcript level
gtf2GRangesList <- function(myfile="/2TB/projects/cleveland/info/Cabili_lincRNAs.gff") {
  gtf <- read.delim(myfile, header=FALSE)
  colnames(gtf) <- c("seqname", "source", "feature", "start", "end", "score", "strand", "frame",      
                     "attributes")
  chronly <- c(1:22, "X", "Y", "MT")
  gtf <- gtf[as.character(gtf$seqname) %in% chronly, ] # Cleanup to remove non-chromosome rows
  
  gene_ids <-  gsub(".*gene_id (.*?);.*", "\\1", gtf$attributes) #get gene_id from attributes column
  transcript_ids<- gsub(".*transcript_id (.*?);.*", "\\1", gtf$attributes) #get transcript_id from attributes column
  
  index<-gene_ids!="" #skip those have no value
  index2<-transcript_ids!=""
  
  gene.gr<-GRanges(seqnames=gtf$seqname[index],
                   ranges=IRanges(gtf$start[index],gtf$end[index]),
                   strand=gtf$strand[index],
                   tx_id=transcript_ids[index],
                   gene_id=gene_ids[index])
  gene.gr.list<-split(gene.gr,gene_ids[index])
  
  transcript.gr<-GRanges(seqnames=gtf$seqname[index2],
                         ranges=IRanges(gtf$start[index2],gtf$end[index2]),
                         strand=gtf$strand[index2],
                         tx_id=transcript_ids[index2],
                         gene_id=gene_ids[index2])
  transcript.gr.list<-split(transcript.gr,transcript_ids[index2])
  
  r<-list()
  r$gene<-gene.gr.list
  r$transcript<-transcript.gr.list
  return(r)
}

lncRNA<-'/scratch/cqs/guoy1/cleveland/lncRNAGFF/Cabili_lincRNAs.gtf'
protein<-'/data/cqs/guoy1/reference/annotation/hg19/Homo_sapiens.GRCh37.63_protein_coding_chr1-22-X-Y-M.gtf'

cat("Loading lncRNA gtf",date())
lnc.gr.list<-gtf2GRangesList(lncRNA)

cat("Loading protein gtf",date())
protein.gr.list<-gtf2GRangesList(protein)

cat("Loading bam file",date())
aln<-readBamGappedAlignments(bamfile)

#
outputcount<-function(tx,aln,outfile){
  cat("Output ",outfile, date())
  counts<-countOverlaps(tx,aln)
  names(counts)<-names(tx)
  write.table(counts,file=outfile,sep="\t")
}

outfile=paste(sample,'lncRNA','gene.txt',sep="_")
outputcount(lnc.gr.list$gene,aln,outfile)

outfile=paste(sample,'lncRNA','isoform.txt',sep="_")
outputcount(lnc.gr.list$transcript,aln,outfile)



outfile=paste(sample,'protein','gene.txt',sep="_")
outputcount(protein.gr.list$gene,aln,outfile)

outfile=paste(sample,'protein','isoform.txt',sep="_")
outputcount(protein.gr.list$transcript,aln,outfile)




setwd('/scratch/cqs/guoy1/cleveland/tophat2_G')





