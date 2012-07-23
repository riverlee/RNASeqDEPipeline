RNASeqDEPipeline (Version 0.1)
================

RNA-Seq differential gene expression pipeline 

## Table of content
  * Usage
  * Required Packages
  * Example
  
## Usage 

```
RNASeqDEPipeline.pl RNASeqDEPipeline.cfg
```

Detailed information about configure file please refer to **RNASeqDEPipeline.cfg** or **example/liver_kind_local.cfg**. A file describing input fastq files has 5 columns, there are **sample**, **group**,**is_paired**,**whichend** and **filename**, please refer to **example/files_local.list** for more information.

## Required Packages
The pipeline is writtein in a combinatin of perl and R. In order to run tophat simutanesouly, package **Parallel::ForkManager** is required. To install it, typing
```
cpan INSTALL Parallel::ForkManager
```
or
```
wget http://search.cpan.org/CPAN/authors/id/D/DL/DLUX/Parallel-ForkManager-0.7.9.tar.gz
tar zxvf Parallel-ForkManager-0.7.9.tar.gz
cd Parallel-ForkManager-0.7.9/
perl Makefile.PL 
make
sudo make install
```
**Installing R packages**
```
#Install packages from CRAN
cran.packages<-c("knitr","markdown","snow","multicore","gplots")
install.packages(cran.packages)

#Install packages from bioconductor
bioconductor.packages<-c("DESeq","baySeq","edgeR","Rsamtools","IRanges","GenomicRanges")
source("http://bioconductor.org/biocLite.R")
biocLite(bioconductor.packages)
```

## Example
We use the RNA-Seq data comes from paper **RNA-seq: An assessment of technical reproducibility and comparison with gene expression arrays. Genome Research 2008 Sep;18(9):1509-17.**. Liver and kidney cell lines with two different concentrations (1.5pm and 3pm) individually were sequenced across the two runs. Here we only used liver and kidney cells with concentration 1.5pm in the second run in our pipleline as an example. The liver samples were placed in lane 1 and lane 7 while the kidney samples were placed in lane 4 and lane 8. The perl codes in the **example/codes** were used to parse the data downloaded from SRA. Here is the steps to use the code
  *download.pl -- download the data from SRA
  *seperate2lane.pl -- Separate reads in a run into each lane. (use fastq-dump first to convert .sra to fastq)
  *change_name.pl -- According to the paper, change the SRRXXXXXX ids into meaningful names
  
  
The configure file and files list description used in this example are **liver_kind_local.cfg** and **files_local.list** in the **example** folder. The output report is named **DEReport.html** containing the running codes and running time in each step (The tophat, cuffdiff, and DE by R are not gitted here).






