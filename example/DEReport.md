Differential expressed gene (DEG) report (RNA-Seq) by RNASeqDEPipeline
========================================================
<div align='right' style='font-weight: bold;'>Created Time: Sun Jul 22 17:21:17 2012</div>
********************************************************
## Table of contents
1. <a href='#M1' target='_self'> Alignment (Tophat)</a>
2. <a href='#M2' target='_self'> DEG by Cuffdiff </a>
3. <a href='#M3' target='_self'> DEG analysis in R </a>
  1. <a href='#M31' target='_self'> Estimate count level expression </a>
  2. <a href='#M32' target='_self'> DEG detection by DESeq, edgeR and baySeq</a>
  3. <a href='#M33' target='_self'> Significant DEG overlaps among different methods </a>
4. <a href='#M4' target='_self'> Session Info </a>

## <a id='M1'></a> 1. Alignment (Tophat) 
**[TopHat](http://tophat.cbcb.umd.edu/)** is a fast splice junction mapper for RNA-Seq reads. It aligns RNA-Seq reads to mammalian-sized genomes using the ultra high-throughput short read aligner **[Bowtie](http://bowtie-bio.sourceforge.net/index.shtml)**, and then analyzes the mapping results to identify splice junctions between exons. 

**Running Time**: 00:05:07:49

**Running Code**:


```r
tophat -p 8  -o tophat/kidneylane4 /home/lij17/accre/gpfs20/data/cqs/guoy1/reference/hg19/bowtie_index/hg19 /2TB/projects/RNA-Seq/normalization/data/rna_normal/lane/kidney1.5pm_run2_lane4.txt
tophat -p 8  -o tophat/kidneylane8 /home/lij17/accre/gpfs20/data/cqs/guoy1/reference/hg19/bowtie_index/hg19 /2TB/projects/RNA-Seq/normalization/data/rna_normal/lane/kidney1.5pm_run2_lane8.txt
tophat -p 8  -o tophat/liverlane1 /home/lij17/accre/gpfs20/data/cqs/guoy1/reference/hg19/bowtie_index/hg19 /2TB/projects/RNA-Seq/normalization/data/rna_normal/lane/liver1.5pm_run2_lane1.txt
tophat -p 8  -o tophat/liverlane7 /home/lij17/accre/gpfs20/data/cqs/guoy1/reference/hg19/bowtie_index/hg19 /2TB/projects/RNA-Seq/normalization/data/rna_normal/lane/liver1.5pm_run2_lane7.txt
```





## <a id='M2'></a> 2. DEG by Cuffdiff 
**[cuffdiff](http://cufflinks.cbcb.umd.edu/manual.html#cuffdiff)** tests for differential expression and is part of **[cufflinks](http://cufflinks.cbcb.umd.edu/)** which assembles transcripts, estimates their abundances, and tests for differential expression and regulation in RNA-Seq samples. 

**Result Files**
  1. [Differentially expressed genes](DE/cuffdiff/gene_exp.diff)
  2. [Differentially expressed transcriprts](DE/cuffdiff/isoform_exp.diff)

**Running Time**: 00:00:04:48

**Running Code**:


```r
cuffdiff --upper-quartile-norm -p 8  -o DE/cuffdiff /home/lij17/accre/gpfs20/data/cqs/guoy1/reference/annotation/hg19/Homo_sapiens.GRCh37.63_protein_coding_chr1-22-X-Y-M.gtf tophat/kidneylane4/accepted_hits.bam,tophat/kidneylane8/accepted_hits.bam  tophat/liverlane1/accepted_hits.bam,tophat/liverlane7/accepted_hits.bam  
```




## <a id='M3'></a> 3. Differential gene expression in R 

### <a id='M31'></a> 3.1 Estimate count level expression 

**Result Files**
  1. [Gene count table](DE/R/gene_expr_count.csv)
  2. [Transcriprt count table](DE/R/transcript_expr_count.csv)

**Running Code**:


```r
# For parallel
suppressPackageStartupMessages(library(multicore))
# Functions defined to calculate count level data for gene/transcript
source("/home/lij17/Dropbox/Documents/RNA-Seq/RNASeqDEPipeline/bam2count.R")
## parameters

## bamfiles are sorted by groups
bamfiles <- c("tophat/kidneylane8/accepted_hits.bam", "tophat/kidneylane4/accepted_hits.bam", 
    "tophat/liverlane1/accepted_hits.bam", "tophat/liverlane7/accepted_hits.bam")
names(bamfiles) <- c("kidneylane8", "kidneylane4", "liverlane1", "liverlane7")
gtffile <- "/home/lij17/accre/gpfs20/data/cqs/guoy1/reference/annotation/hg19/Homo_sapiens.GRCh37.63_protein_coding_chr1-22-X-Y-M.gtf"
threads <- 2
logfile <- "run.log"
geneoutput <- "DE/R/gene_expr_count.csv"
transcriptoutput <- "DE/R/transcript_expr_count.csv"

# Step 1) load gtf
write.log(msg = "Load GTF file...", file = logfile)
system.time(gr <- gtf2GRangesList(myfile = gtffile))
```

```
##    user  system elapsed 
## 188.836   0.964 190.185 
```

```r

# Step 2) Loop bamfiles, get count level for both gene and transcript
write.log(msg = "Start to get count level expression from bamfiles", file = logfile)
system.time(count.list <- mclapply(bamfiles, function(f) {
    write.log(msg = paste("Getting count levele expression from file '", f, 
        "'", sep = ""), file = logfile)
    aln = readBamGappedAlignments(f)
    r <- mclapply(gr, function(x) {
        getCounts(aln, grangelist = x)
    }, mc.cores = 2)
}, mc.cores = threads))
```

```
##    user  system elapsed 
##  90.085   7.536  81.243 
```

```r

gene.count <- sapply(count.list, function(x) x$gene)
transcript.count <- sapply(count.list, function(x) x$transcript)

# Step 3) Write out result
write.log(msg = "Write out count level data into files", file = logfile)
write.csv(gene.count, file = geneoutput)
write.csv(transcript.count, file = transcriptoutput)
```




###  3.2 Differential gene expression by DESeq, edgeR, baySeq and TSPM <a id='M32'></a>

**Method brief description**
  1. **DESeq** --  Simon Anders, Wolfgang Huber: Differential expression analysis for sequence count data. Genome Biology 11 (2010) R106, [doi:10.1186/gb-2010-11-10-r106](http://dx.doi.org/10.1186/gb-2010-11-10-r106)/ [R package at bioconductor](http://bioconductor.org/packages/release/bioc/html/DESeq.html)
  2. **edgeR** -- Robinson MD, McCarthy DJ, Smyth GK: edgeR: a Bioconductor package for differential expression analysis of digital gene expression data. Bioinformatics (2010), [doi:10.1093/bioinformatics/btp616](http://dx.crossref.org/10.1093%2Fbioinformatics%2Fbtp616) /  [R package at bioconductor](http://bioconductor.org/packages/release/bioc/html/edgeR.html)
  3. **baySeq** -- Thomas J Hardcastle, Krystyna A Kelly: baySeq: Empirical Bayesian methods for identifying differential expression in sequence count data. BMC Bioinformatics (2010), [doi:10.1186/1471-2105-11-422](http://dx.doi.org/10.1186/1471-2105-11-422) / [R package at bioconductor](http://bioconductor.org/packages/release/bioc/html/baySeq.html)
  4. **TSPM** -- Paul L. Auer, Rebecca W. Doerge: A Two-Stage Poisson Model for Testing RNA-Seq Data. Statistical Applications in Genetics and Molecular Biology (2011), [doi: 10.2202/1544-6115.1627](http://dx.doi.org/10.2202/1544-6115.1627) / [R Code](http://www.stat.purdue.edu/~doerge/software/TSPM.R)

**Result Files**
  1. [DE genes by R](DE/R/de_gene_expr_count.csv)
  2. [DE transcriprts by R](DE/R/de_transcript_expr_count.csv)

**Running Code**:



```r
source("/home/lij17/Dropbox/Documents/RNA-Seq/RNASeqDEPipeline/DEUtil.R")
group <- c("1", "1", "2", "2")
mc.cores <- 4
degeneoutfile <- "DE/R/de_gene_expr_count.csv"
detranscriptoutfile <- "DE/R/de_transcript_expr_count.csv"

## Do differential gene expression
countData <- read.csv("DE/R/gene_expr_count.csv", row.names = 1)
system.time(wrapFourDE(countData, group, degeneoutfile, mc.cores = mc.cores))
```

```
##    user  system elapsed 
## 121.620   9.369 319.182 
```

```r

## Do differentail transcript expression
countData <- read.csv("DE/R/transcript_expr_count.csv", row.names = 1)
system.time(wrapFourDE(countData, group, detranscriptoutfile, mc.cores = mc.cores))
```

```
##     user   system  elapsed 
##  471.789    6.476 1980.763 
```





### 3.3 Significant DEG overlaps among different methods <a id='M33'></a>
This module lists all the significant genes and transcript among all the methods at FDR<=0.1
**Result Files**
  1. [Overlapped DE genes among all the methods](DE/R/overlapped_de_gene_expr_count.csv)
  2. [Overlapped DE transcriprts among all the methods](DE/R/overlapped_de_transcript_expr_count.csv)
  3. [All DE genes](DE/R/all_de_gene_expr_count.csv)
  4. [All DE transcriprts](DE/R/all_de_transcript_expr_count.csv)
  
**Running Code**:


```r
source("/home/lij17/Dropbox/Documents/RNA-Seq/RNASeqDEPipeline/seeOverlaps.R")
# Overlapped DE genes
cuffdiffFile.gene <- "DE/cuffdiff/gene_exp.diff"
rdeFile.gene <- "DE/R/de_gene_expr_count.csv"
allDE.gene.output <- "DE/R/all_de_gene_expr_count.csv"
overlapDE.gene.output <- "DE/R/overlapped_de_gene_expr_count.csv"
system.time(seeOverlap(cuffdiffFile = cuffdiffFile.gene, rdeFile = rdeFile.gene, 
    output.all = allDE.gene.output, output.overlap = overlapDE.gene.output, 
    type = "genes"))
```

```
## There are 2913 overlapped genes among all the methods at FDR<=0.1
```

![plot of chunk overlap](figure/overlap1.png) ![plot of chunk overlap](figure/overlap2.png) ![plot of chunk overlap](figure/overlap3.png) 

```
##    user  system elapsed 
##   3.564   0.056   3.628 
```





```r
# Overlapped DE transcript
cuffdiffFile.transcript <- "DE/cuffdiff/isoform_exp.diff"
rdeFile.transcript <- "DE/R/de_transcript_expr_count.csv"
allDE.transcript.output <- "DE/R/all_de_transcript_expr_count.csv"
overlapDE.transcript.output <- "DE/R/overlapped_de_transcript_expr_count.csv"
system.time(seeOverlap(cuffdiffFile = cuffdiffFile.transcript, rdeFile = rdeFile.transcript, 
    output.all = allDE.transcript.output, output.overlap = overlapDE.transcript.output, 
    type = "transcripts"))
```

```
## There are 4333 overlapped transcripts among all the methods at FDR<=0.1
```

![plot of chunk overlap3](figure/overlap31.png) ![plot of chunk overlap3](figure/overlap32.png) ![plot of chunk overlap3](figure/overlap33.png) 

```
##    user  system elapsed 
##  28.678   0.124  28.854 
```



## 4. Session Info <a id='M4'></a>


```r
sessionInfo()
```

```
## R version 2.14.1 (2011-12-22)
## Platform: x86_64-pc-linux-gnu (64-bit)
## 
## locale:
##  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
##  [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
##  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
##  [7] LC_PAPER=C                 LC_NAME=C                 
##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
## [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
## 
## attached base packages:
## [1] grid      stats     graphics  grDevices utils     datasets  methods  
## [8] base     
## 
## other attached packages:
##  [1] gplots_2.11.0       MASS_7.3-19         KernSmooth_2.23-8  
##  [4] caTools_1.13        bitops_1.0-4.1      gdata_2.11.0       
##  [7] gtools_2.7.0        snow_0.3-10         baySeq_1.8.3       
## [10] edgeR_2.4.6         limma_3.10.3        DESeq_1.6.1        
## [13] locfit_1.5-8        Biobase_2.14.0      Rsamtools_1.6.3    
## [16] Biostrings_2.22.0   GenomicRanges_1.6.7 IRanges_1.12.6     
## [19] multicore_0.1-7     markdown_0.5.2      knitr_0.6.3        
## 
## loaded via a namespace (and not attached):
##  [1] annotate_1.32.3       AnnotationDbi_1.16.19 BSgenome_1.22.0      
##  [4] DBI_0.2-5             digest_0.5.2          evaluate_0.4.2       
##  [7] formatR_0.6           genefilter_1.36.0     geneplotter_1.32.1   
## [10] lattice_0.20-6        plyr_1.7.1            RColorBrewer_1.0-5   
## [13] RCurl_1.6-6           RSQLite_0.11.1        rtracklayer_1.14.4   
## [16] splines_2.14.1        stringr_0.6           survival_2.36-14     
## [19] tools_2.14.1          XML_3.9-4             xtable_1.7-0         
## [22] zlibbioc_1.0.1       
```





<div style="background-color: #EEE;font-weight: bold;" align='center'>Produced by <a href="https://github.com/riverlee/RNASeqDEPipeline">RNASeqDEPipeline</a> (Jiang Li)</div>
</body>

