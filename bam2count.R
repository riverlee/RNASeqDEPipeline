suppressPackageStartupMessages(library(Rsamtools))
suppressPackageStartupMessages(library(IRanges))
suppressPackageStartupMessages(library(GenomicRanges))

#################################
# Start function defination
# Fun1: Return GRangesList object in both gene level and transcript level
gtf2GRangesList <- function(myfile=null) {
  gtf <- read.delim(myfile, header=FALSE)
  colnames(gtf) <- c("seqname", "source", "feature", "start", "end", "score", "strand", "frame",      
                     "attributes")
  #chronly <- c(1:22, "X", "Y", "MT")
  #gtf <- gtf[as.character(gtf$seqname) %in% chronly, ] # Cleanup to remove non-chromosome rows
  
  gene.ids <-  gsub(".*gene_id (.*?);.*", "\\1", gtf$attributes) #get gene_id from attributes column
  #gene.name <- gsub(".*gene_name (.*?);.*", "\\1", gtf$attributes) #get gene_name from attributes column
  transcript.ids<- gsub(".*transcript_id (.*?);.*", "\\1", gtf$attributes) #get transcript_id from attributes column
  
  gene.index<-gene.ids!=""                    #skip those have no value
  transcript.index<-transcript.ids!=""   #skip those have no value
  
  gene.gr<-GRanges(seqnames=gtf$seqname[gene.index],
                   ranges=IRanges(gtf$start[gene.index],gtf$end[gene.index]),
                   strand=gtf$strand[gene.index],
                   tx_id=transcript.ids[gene.index],
                   gene_id=gene.ids[gene.index])
  gene.gr.list<-split(gene.gr,gene.ids[gene.index])
  
  transcript.gr<-GRanges(seqnames=gtf$seqname[transcript.index],
                         ranges=IRanges(gtf$start[transcript.index],gtf$end[transcript.index]),
                         strand=gtf$strand[transcript.index],
                         tx_id=transcript.ids[transcript.index],
                         gene_id=gene.ids[transcript.index])
  transcript.gr.list<-split(transcript.gr,transcript.ids[transcript.index])
  
  r<-list()
  r$gene<-gene.gr.list
  r$transcript<-transcript.gr.list
  return(r)
}

# Fun2: write out message to log file
write.log<-function(msg=NULL,file=NULL,append=TRUE){
	cat("[Info_R][",date(),"]:",msg,"\n",sep="",file=file,append=append)
}

# Fun3: get count based on bam oject  and grangelist
getCounts<-function(aln,grangelist){
  counts<-countOverlaps(grangelist,aln)
  names(counts)<-names(grangelist)
  return(counts)
}
