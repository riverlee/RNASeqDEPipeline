#!/usr/bin/Rscript
#############################################
#Author: Jiang Li
#email: riverlee2008@gmail.com
#Creat Time: 
#Vanderbilt Center for Quantitative Sciences
#############################################
suppressPackageStartupMessages(library(gplots))
#suppressPackageStartupMessages(library("googleVis"))

print.table<-function(data=NULL,width=1000,height=520,pageSize=50,...){
  PopTable <- gvisTable(data, options = list(width = width, height = height,pageSize=pageSize,page = "enable"))
  print(PopTable,'chart')
}


seeOverlap<-function(cuffdiffFile=NULL,rdeFile=NULL,output.all=NULL,output.overlap=NULL,type='Genes',printTable=TRUE){
  cuffdiff.de<-read.delim(cuffdiffFile)
  cuffdiff.de.sub<-cuffdiff.de[,c(3,4,10,12,13)]
  colnames(cuffdiff.de.sub)<-c('Gene','Locus','cuffdiff_log2FC','cuffdiff_pvalue','cuffdiff_FDR')
  rownames(cuffdiff.de.sub)<-cuffdiff.de$test_id
  r.de<-read.csv(rdeFile,row.names=1)
  overlapped.ids<-intersect(cuffdiff.de$test_id,r.de$ID)
  cuffdiff.de.sub<-cuffdiff.de.sub[rownames(cuffdiff.de.sub) %in% overlapped.ids,]
  r.de<-r.de[rownames(r.de) %in% overlapped.ids,2:ncol(r.de)]
  
  #Construct new data set
  result.data.frame<-data.frame(ID=overlapped.ids[order(overlapped.ids)],cuffdiff.de.sub[order(rownames(cuffdiff.de.sub)),],r.de[order(rownames(r.de)),])
  rownames(result.data.frame)<-result.data.frame$ID
  
  # write out all merged result
  write.csv(result.data.frame,output.all)
  
  # exact fdrs
  index <- grepl("_FDR$", colnames(result.data.frame), perl = TRUE)
  fdrs <- result.data.frame[, index]
  fdrs[is.na(fdrs)] <- 1
  colnames(fdrs) <- gsub("_FDR","",colnames(fdrs))
  
  # Output overlapped genes at FDR <=0.1
  sig.list<-apply(fdrs,2,function(x) rownames(fdrs)[x<=0.1])
  overlapped.genes<-character()
  if(length(sig.list)>=3){
  	overlapped.genes<-intersect(sig.list[[1]],sig.list[[2]])
  	for(i in 3:length(sig.list)){
  		overlapped.genes<-intersect(overlapped.genes,sig.list[[i]])
  	}
  }else if(length(sig.list)==2){
  	overlapped.genes<-intersect(sig.list[[1]],sig.list[[2]])	
  }else{
  	overlapped.genes<-sig.list[[1]]
  }
  overlapped.result<-result.data.frame[overlapped.genes,]
  cat("There are ",length(overlapped.genes)," overlapped ",type," among all the methods at FDR<=0.1",sep="")
  write.csv(result.data.frame[overlapped.genes,],output.overlap)
  #print.table(result.data.frame[1:100,])
  
  # plot venn diagram
  fdrcutoffs<-c(0.01,0.05,0.1)
  for (cutoff in fdrcutoffs){   
  	sig.list<-apply(fdrs,2,function(x) rownames(fdrs)[x<=cutoff])
  	venn(sig.list)
  	title(paste("Number of significant ",type," at FDR<=",cutoff,sep=""))
  }  
}





