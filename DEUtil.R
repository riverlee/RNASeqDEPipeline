suppressPackageStartupMessages(library(DESeq))
suppressPackageStartupMessages(library(edgeR))
suppressPackageStartupMessages(library(baySeq))
suppressPackageStartupMessages(library(snow))
suppressPackageStartupMessages(library(multicore))


##############################################
# Fuctions wrap four methods together
# will be invoked by wrapFourDE

callFun <- function(fun, countData, group) {
    tryCatch(eval(fun(countData, group)),
        error=function(e) return("error"))
}

wrapFourDE <- function(countData,group, outfile,mc.cores=2) {
    funs <- c(run.DESeq,run.baySeq, run.edgeR,run.TSPM)
    result <- mclapply(funs, callFun, countData = countData, group = group, 
        mc.cores = mc.cores)
    names(result) <- c("DESeq", "baySeq","edgeR", "TSPM")
    
    error.index<-unlist(lapply(result,function(x){
      ifelse(length(x)==1 && x == 'error',TRUE,FALSE)
    }))
    
    if(sum(error.index)){
      if(sum(error.index)!=length(funs)){
        cat("[Warning]: Test of '",paste(names(result)[error.index], collapse="', '"),"' failed",sep="")
        result<-result[!error.index]
      }else{
        stop("All the tests are failed")
      }
    }
    
    result.data.frame <- data.frame(ID = rownames(result[[1]])[order(rownames(result[[1]]))], stringsAsFactors = FALSE)
  
    for (d in result) {
        result.data.frame <- cbind(result.data.frame, d[order(rownames(d)), 
            ])
    }

    colnames(result.data.frame)[2:ncol(result.data.frame)] <- paste(rep(names(result), 
        each = 3), colnames(result.data.frame)[2:ncol(result.data.frame)], sep = "_")
    write.csv(result.data.frame, outfile)
}

##############################################
# Functions for DE of RNA-Seq

# Package: DESeq

run.DESeq<-function(countData, group, locfunc=median,method="per-condition",fitType='local',sharingMode='maximum',...){
  #Check 
  conditions<-unique(group)
  if(length(conditions)!=2){
  	stop("'group' doesn't contain two conditions")
  }
  if(ncol(countData) !=length(group)){
  	stop("'group' label doesn't meet with 'countData' columns")
  }
  
  #main part
  y_DESeq<-newCountDataSet( countData, group )
  y_DESeq<-estimateSizeFactors( y_DESeq,locfunc=locfunc )
  y_DESeq<-estimateDispersions( y_DESeq , method=method, fitType=fitType,sharingMode)
  #y_DESeq<-estimateVarianceFunctions(y_DESeq)
  DESeq.r<-nbinomTest( y_DESeq, conditions[1],conditions[2])

  final.r<-data.frame(log2FC=DESeq.r$log2FoldChange, pvalue=DESeq.r$pval,FDR=DESeq.r$padj)
  rownames(final.r)<-DESeq.r$id
  final.r<-final.r[order(final.r$FDR,final.r$pvalue,-abs(final.r$log2FC)),]
  final.r  
}


# Package: edgeR
# method could be TMM, upperquartile, RLE
run.edgeR<-function(countData,group,method="TMM",quantile=0.75,...){
  #Check 
  conditions<-unique(group)
  if(length(conditions)!=2){
  	stop("'group' doesn't contain two conditions")
  }
  if(ncol(countData) !=length(group)){
  	stop("'group' label doesn't meet with 'countData' columns")
  }

  #main part	
  y<-DGEList(countData, group=group )
  n.top<-nrow(countData) #return all the genes
  edgeR.r<-topTags( exactTest( 
             estimateTagwiseDisp( 
               estimateCommonDisp( 
                 calcNormFactors(y, method="TMM",p=quantile) 
               )
             ) 
           ), 
           n.top)$table
  
  final.r<-data.frame(log2FC=edgeR.r$logFC, pvalue=edgeR.r$PValue,FDR=edgeR.r$FDR)
  rownames(final.r)<-rownames(edgeR.r)
  final.r<-final.r[order(final.r$FDR,final.r$pvalue,-abs(final.r$log2FC)),]
  final.r                                       
}


# Pacakge: baySeq
# estimationType could be one of quantile, total, or edgeR (TMM: Trimmed Mean of M-vales)
# estimation could be one of QL, ML and edgeR
run.baySeq<-function(countData,group,estimationType="quantile",quantile=0.75,estimation="ML",samplesize=1000,parallel=TRUE,threads=4,...){
  #Check 
  conditions<-unique(group)
  if(length(conditions)!=2){
  	stop("'group' doesn't contain two conditions")
  }
  if(ncol(countData) !=length(group)){
  	stop("'group' label doesn't meet with 'countData' columns")
  }
  
  #main part
  cl<-NULL
  if(parallel){
    # IF FAILED, will set cl=null (no parallele)
    tryCatch(cl<-makeCluster(threads,"SOCK"),error=function(e) return(cl<<-NULL))
  }
  countData<-as.matrix(countData)
  n.top<-nrow(countData) #return all the genes
  conds <- list(NDE = rep(1, length(group)),DE = group)
  CD<-new("countData", data = countData, replicates = group, groups = conds)
  CD@libsizes <- getLibsizes(CD, estimationType=estimationType,quantile=quantile)
  CD <- getPriors.NB(CD, samplesize = samplesize, estimation = estimation, equalDispersions=TRUE, cl = cl)
  CD <- getLikelihoods.NB(CD, pET = 'BIC', cl = cl)
  baySeq.r<-topCounts(CD, group = "DE", number=n.top,normaliseData=TRUE)
  if(!is.null(cl)){
  	stopCluster(cl)
  }
  # normalisation then get log2foldchange
  samplesA<-group==conditions[1]
  samplesB<-group==conditions[2]
  Adata <- CD@data[, samplesA]
  Bdata <- CD@data[, samplesB]
  Adata <- Adata/CD@libsizes[samplesA] * mean(CD@libsizes)
  Bdata <- Bdata/CD@libsizes[samplesB] * mean(CD@libsizes)
  logFC<-log2(rowMeans(Bdata)/rowMeans(Adata))
  
  final.r<-data.frame(log2FC=logFC[rownames(baySeq.r)], Likelihood=baySeq.r$Likelihood,FDR=baySeq.r$FDR)
  rownames(final.r)<-rownames(baySeq.r)
  final.r<-final.r[order(final.r$FDR,-final.r$Likelihood,-abs(final.r$log2FC)),]
  final.r
}


# Pacakge: TSPM source code from http://www.stat.purdue.edu/~doerge/software/TSPM.R
# estimationType could be one of quantile, total, or edgeR (TMM: Trimmed Mean of M-vales)
# estimation could be one of QL, ML and edgeR
run.TSPM<-function(countData,group,estimationType="quantile",quantile=0.75,parallel=TRUE,threads=4,...){
  #Check 
  conditions<-unique(group)
  if(length(conditions)!=2){
  	stop("'group' doesn't contain two conditions")
  }
  if(ncol(countData) !=length(group)){
  	stop("'group' label doesn't meet with 'countData' columns")
  }
  
  #main part
  cl<-NULL
  if(parallel){
    tryCatch(cl<-makeCluster(threads,"SOCK"),error=function(e) return(cl<<-NULL))
  }
  countData<-as.matrix(countData)
  n.top<-nrow(countData) #return all the genes
  conds <- list(NDE = rep(1, length(group)),DE = group)
  CD<-new("countData", data = countData, replicates = group, groups = conds)
  lib.size_UpQ<-getLibsizes(CD, estimationType=estimationType,quantile=quantile)
  TSPM.r<-TSPM(countData, x1=group, x0=conds$NDE, lib.size_UpQ,cl=cl)
  if(!is.null(cl)){
  	stopCluster(cl)
  }
  
  final.r<-data.frame(log2FC=TSPM.r$result$log.fold.change,pvalue=TSPM.r$result$pvalues,FDR=TSPM.r$result$padj)
  rownames(final.r)<-rownames(countData)
  final.r<-final.r[order(final.r$FDR,final.r$pvalue,-abs(final.r$log2FC)),]
  final.r 
}


## From: http://www.stat.purdue.edu/~doerge/software/TSPM.R
##-------------------------------------------------------------------
## 	Name: TSPM.R
##  	R code for the paper by Paul L. Auer and R.W. Doerge:
## 	"A Two-Stage Poisson Model for Testing RNA-Seq Data"
## 	Date: February 2011
## 	Contact: Paul Auer 		plivermo@fhcrc.org
##		 R.W. Doerge 		doerge@purdue.edu


## 	Example:
## 	counts <- matrix(0, nrow=1000, ncol=10)
## 	for(i in 1:1000){
##		lambda <- rpois(n=1, lambda=10)
##		counts[i,] <- rpois(n=10, lambda=lambda)
## 	}
## 	x1 <- gl(n=2, k=5, labels=c("T", "C"))
## 	x0 <- rep(1, times=10)
## 	lib.size <- apply(counts,2,sum)
##	result <- TSPM(counts, x1, x0, lib.size)
##---------------------------------------------------------------------


#######################################################################
###### The TSPM function ##############################################
#######################################################################


TSPM <- function(counts, x1, x0, lib.size, alpha.wh=0.05,cl=NULL){

## Input:
#counts: 	a matrix of RNA-Seq gene counts (genes are rows, samples are columns)
#x1: 		a vector of treatment group factors (under the alternative hypothesis)
#x0: 		a vector of treatment group factors (under the null hypothesis)
#lib.size: 	a vector of RNA-Seq library sizes. This could simply be obtained
#          	by specifying lib.size <- apply(counts,2,sum). It may also be any other
#          	appropriate scaling factor.
#alpha.wh:	the significance threshold to use for deciding whether a gene is overdispersed.
#               Defaults to 0.05.


## Output:
#log.fold.change:		a vector containing the estimated log fold changes for each gene
#pvalues: 			a vector containing the raw p-values testing differential expression for each gene.
#index.over.disp: 	a vector of integer values containing the indices of the over-dispersed genes.
#index.not.over.disp:	a vector of integer values containing the indices of the non-over-dispersed genes.
#padj:			a vector containing the p-values after adjusting for multiple testing using the 
#				method of Benjamini-Hochberg


######## The main loop that fits the GLMs to each gene ########################

### Initializing model parameters ####
n <- dim(counts)[1]
per.gene.disp <- NULL
LRT <- NULL
score.test <- NULL
LFC <- NULL

###### Fitting the GLMs for each gene #################
# Old one
#	for(i in 1:n){
#		### Fit full and reduced models ###
#		model.1 <- glm(as.numeric(counts[i,]) ~ x1, offset=log(lib.size), family=poisson)
#		model.0 <- glm(as.numeric(counts[i,]) ~ x0, offset=log(lib.size), family=poisson)

#		### Obtain diagonals of Hat matrix from the full model fit ###
#		hats <- hatvalues(model.1)

#		### Obtain Pearson overdispersion estimate ####
#		per.gene.disp[i] <- sum(residuals(model.1, type="pearson")^2)/model.1$df.residual

#		### Obtain Likelihood ratio statistic ####
#		LRT[i] <- deviance(model.0)-deviance(model.1)

#		### Obtain score test statistic ####
#		score.test[i] <- 1/(2*length(counts[i,])) * sum(residuals(model.1, type="pearson")^2 - ((counts[i,] - hats*model.1$fitted.values)/model.1$fitted.values))^2
#		
#		### Obtain the estimated log fold change ###
#		LFC[i] <- -model.1$coef[2]
#	}
	
	#updated one by Jiang Li: Wed 18 Jul 2012 11:16:26 PM CDT 
	fill.value<-function(i,x1,x0,counts,lib.size){
		### Fit full and reduced models ###
		model.1 <- glm(as.numeric(counts[i,]) ~ x1, offset=log(lib.size), family=poisson)
		model.0 <- glm(as.numeric(counts[i,]) ~ x0, offset=log(lib.size), family=poisson)

		### Obtain diagonals of Hat matrix from the full model fit ###
		hats <- hatvalues(model.1)

		### Obtain Pearson overdispersion estimate ####
		per.gene.disp <- sum(residuals(model.1, type="pearson")^2)/model.1$df.residual

		### Obtain Likelihood ratio statistic ####
		LRT<- deviance(model.0)-deviance(model.1)

		### Obtain score test statistic ####
		score.test<- 1/(2*length(counts[i,])) * sum(residuals(model.1, type="pearson")^2 - ((counts[i,] - hats*model.1$fitted.values)/model.1$fitted.values))^2
	
		### Obtain the estimated log fold change ###
		LFC<- -model.1$coef[2]
		t<-c(per.gene.disp,LRT,score.test,LFC)
		#names(t)<-c("per.gene.disp","LRT","score.test","LFC")
		t
	}
	
	if(is.null(cl)){
		tmp<-lapply(1:n,fill.value,x1=x1,x0=x0,counts=counts,lib.size=lib.size)
	}else{
		tmp<-parLapply(cl,1:n,fill.value,x1=x1,x0=x0,counts=counts,lib.size=lib.size) #use parallele
	}
	
	per.gene.disp<-unlist(lapply(tmp,function(x) x[1]))
	LRT <- unlist(lapply(tmp,function(x) x[2]))
	score.test <- unlist(lapply(tmp,function(x) x[3]))
	LFC <- unlist(lapply(tmp,function(x) x[4]))

## Initialize parameters for Working-Hotelling bands around the score TSs ###
qchi <- qchisq(df=1, (1:n-0.5)/n)
MSE <- 2
UL <- NULL

#### Obtain the upper boundary of the WH bands #######################################
xbar <- mean(qchi)
bottom <- sum((qchi-xbar)^2)
top <- (qchi-xbar)^2
s <- sqrt(MSE*(1/n) + (top/bottom))
W <- sqrt(2*qf(df1=1, df2=n-1, p=1-(alpha.wh/n)))
UL <- pmax(qchi + W*s,1)

###### Obtain the indices of the over-dispersed and not-over-dispersed genes, respectively ##########

cutoff <- min(which(sort(score.test)-UL > 0))
temp <- cutoff-1 + seq(cutoff:length(score.test))
over.disp <- which(score.test %in% sort(score.test)[temp])
not.over.disp <- setdiff(1:length(score.test), over.disp)

###### Compute p-values ####################################
model.1 <- glm(as.numeric(counts[1,]) ~ x1, offset=log(lib.size), family=poisson) #added by Jiang Li
p.f <- pf(LRT[over.disp]/per.gene.disp[over.disp], df1=1, df2=model.1$df.residual, lower.tail=FALSE)
p.chi <- pchisq(LRT[not.over.disp], df=1, lower.tail=FALSE)
p <- NULL
p[over.disp] <- p.f
p[not.over.disp] <- p.chi

##### Adjust the p-values using the B-H method ####################
p.bh.f <- p.adjust(p.f, method="BH")
p.bh.chi <- p.adjust(p.chi, method="BH")
final.p.bh.tagwise <- NULL
final.p.bh.tagwise[over.disp] <- p.bh.f
final.p.bh.tagwise[not.over.disp] <- p.bh.chi

### Output ###
#list(log.fold.change=LFC, pvalues=p, index.over.disp=over.disp, index.not.over.disp=not.over.disp,
# padj=final.p.bh.tagwise)
return( list(result=data.frame(log.fold.change=LFC, pvalues=p, padj=final.p.bh.tagwise), 
        index.over.disp=over.disp, index.not.over.disp=not.over.disp)  )
}

