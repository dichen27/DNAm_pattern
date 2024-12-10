rm(list=ls())
library(minfi)
setwd("/Users/dichen/Documents/PPMI/preprocess/")
datadir <- "/Users/dichen/Documents/PPMI/ppmi_120_methylation_profiling_idat_20180312/Project120_IDATS_n524final_toLONI_030718/"

#setwd("/Users/dichen/Documents/PPMI_Running/PPMI/preprocess/")
#datadir <- "/Users/dichen/Documents/PPMI_Running/PPMI/ppmi_120_methylation_profiling_idat_20180312/Project120_IDATS_n524final_toLONI_030718/"

## read data
p1<-read.metharray.sheet( "/Users/dichen/Documents/PPMI_Running/PPMI/ppmi_120_methylation_profiling_idat_20180312/Project120_IDATS_n524final_toLONI_030718/PPMI_information.csv", recursive=T)



pd <- p1
RGset <- read.metharray.exp(base =datadir,targets=pd, verbose=TRUE,force=T)
#RGset <- read.metharray.exp(targets=pd)

#read.metharray.exp(base = NULL, targets = pd, extended = FALSE, recursive = FALSE, verbose = FALSE, force = FALSE)

save(RGset,file="./RGset_PPMI.rda") 

##### test example

if(require(minfiData)) {
  baseDir <- system.file("extdata", package = "minfiData")
  RGset <- read.metharray.exp(file.path(baseDir, "5723646052"))
}

if(require(minfiData)) {
  baseDir <- system.file("extdata", package = "minfiData")
  RGset <- read.metharray.exp(file.path(baseDir, "5723646052"))
}





####  IMAGEN code
fileName <- "./rdas/RGset.rda"
if (file.exists(fileName)) {
    load(fileName)
} else {
	############
	## wave 1 ##
	############
	datadir <- "/home/imagen/BL Illumina 450k Wave1"
	p1=read.metharray.sheet(datadir, "IMGN1_METH1_PL21_1x_121228.v2.csv", recursive=T)
	p1=rbind(p1, read.metharray.sheet(datadir, "IMGN1_PL20_1x_121228.v2.csv", recursive=T))
	p1=rbind(p1, read.metharray.sheet(datadir, "IMGN1_PL2[2345]_1x_130114.v2.csv$", recursive=T))
	p1=rbind(p1, read.metharray.sheet(datadir, "IMGN1_PL27_1x_130121.v2.csv", recursive=T))
	p1=rbind(p1, read.metharray.sheet(datadir, "IMGN1_SVK_PL26_1x_130121.v2.csv", recursive=T))
	p1=rbind(p1, read.metharray.sheet(datadir, "METH1_IMGN1_SVK_PL3_1x_121206.v2.csv", recursive=T))
	p1 <- p1[!(p1[,"Sample.ID"]==""),]

	############
	## wave 2 ##
	############
    datadir <- "/home/imagen/BL Illumina 450k Wave2"
    p2=read.metharray.sheet(datadir, "wave2-0000-0200.v2.csv", recursive=T)
    p2=rbind(p2, read.metharray.sheet(datadir, "wave2-0201-0400.v2.csv", recursive=T))
    p2=rbind(p2, read.metharray.sheet(datadir, "wave2-0401-0600.v2.csv", recursive=T))
    p2=rbind(p2, read.metharray.sheet(datadir, "wave2-0601-0800.v2.csv", recursive=T))
    p2=rbind(p2, read.metharray.sheet(datadir, "wave2-0801-1000.v2.csv", recursive=T))
    p2=rbind(p2, read.metharray.sheet(datadir, "wave2-1001-1178.v2.csv", recursive=T))
    p2$Notes <- ""

	p1  <- p1[,match(colnames(p2), colnames(p1))]
	stopifnot(identical(colnames(p1),colnames(p2)))
	plates <- rbind(p1,p2)
	Index <- plates[,"Notes"]=="Methylation control"
	plates <- plates[!Index,]
	Index <- plates[,"Sample.ID"]=="IMGN1_901.1_25R04C02" #gender discrepency
	plates <- plates[!Index,]
	plates$platform <- "450k"
	plates <- plates[!is.na(plates$IMAGEN_PSC1),]
	
	############
	## 850K   ##
	############	
    datadir <- "/home/imagen/FU2 Illumina 850k Wave1"
    p3=read.metharray.sheet(datadir, "850K_wave1.csv", recursive=T)
    p3 <- p3[!is.na(p3[,1]),]

    datadir <- "/home/imagen/FU2 Illumina 850k Wave 2"
    p4=read.metharray.sheet(datadir, "850K_wave2.csv", recursive=T)
    p4 <- p4[!is.na(p4[,1]),]
	p4 <- p4[,match(colnames(p3), colnames(p4))]
	p3 <- rbind(p3,p4)
	p3$Notes <- ""
	p3$platform <- "850k"

	plates <- plates[,colnames(plates) %in% colnames(p3)]
	p3 <- p3[,colnames(p3) %in% colnames(plates)]
	p3 <- p3[,match(colnames(plates),colnames(p3))]
	stopifnot(identical(colnames(p3),colnames(plates)))

	pd <- plates
	RGset = read.metharray.exp(targets=pd, verbose=TRUE)
	save(RGset,file="./rdas/RGset_450k.rda")
	RGset1 <- RGset
	
	pd <- p3
	RGset = read.metharray.exp(targets=pd, verbose=TRUE,force=T)
	save(RGset,file="./rdas/RGset_850k.rda")

	rgSet <- combineArrays(RGset1, RGset, outType="IlluminaHumanMethylationEPIC")
	RGset <- rgSet
	save(RGset,file=fileName)
}



