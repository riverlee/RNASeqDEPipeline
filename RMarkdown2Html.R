#!/usr/bin/Rscript
#############################################
#Author: Jiang Li
#email: riverlee2008@gmail.com
#Creat Time: Wed 04 Jul 2012 06:31:12 PM CDT 
#Vanderbilt Center for Quantitative Sciences
#############################################

library(knitr)
library(markdown)

args<-commandArgs(trailingOnly=TRUE)
knit2html(args[1])


