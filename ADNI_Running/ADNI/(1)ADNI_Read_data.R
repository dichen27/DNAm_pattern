rm(list=ls())
library(minfi)

## mac
#setwd("/Users/dichen/Documents/ADNI/")
#datadir <- "/Users/dichen/Documents/ADNI/ADNI_iDAT_files/"

## FD-server
setwd("/share/home1/chendi/Documents/ADNI/")
datadir <- "/share/home1/chendi/Documents/ADNI/ADNI_iDAT_files/"

## read data
#pd<-read.metharray.sheet( "./ADNI_information.csv", recursive=T)

RGset <- read.metharray.exp(base =datadir,targets=NULL, verbose=TRUE,force=T)
#RGset <- read.metharray.exp(targets=pd)

#read.metharray.exp(base = NULL, targets = pd, extended = FALSE, recursive = FALSE, verbose = FALSE, force = FALSE)

save(RGset,file="./RGset_ADNI_1919.rda") 

##### test example

if(require(minfiData)) {
  baseDir <- system.file("extdata", package = "minfiData")
  RGset <- read.metharray.exp(file.path(baseDir, "5723646052"))
}

if(require(minfiData)) {
  baseDir <- system.file("extdata", package = "minfiData")
  RGset <- read.metharray.exp(file.path(baseDir, "5723646052"))
}


