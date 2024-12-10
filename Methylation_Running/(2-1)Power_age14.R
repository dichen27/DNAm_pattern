rm(list=ls())
setwd('/Users/dichen/Documents/Methylation_Running/')

#install.packages("R.matlab")
library(R.matlab)
num <- 372582

## read data
t14_source<-readMat('data_regressed/t14_full_regress.mat')
t14<-t14_source$t14.full.regress
t14[1:3,1:3]
Raw<-t14[,-1]
Raw[1:3,1:3]
nrow(Raw)
ncol(Raw)

#install.packages("BiocManager")
#BiocManager::install("WGCNA")

library(WGCNA)
allowWGCNAThreads()
options(stringsAsFactors=FALSE)
gsg <- goodSamplesGenes(as.matrix(Raw), verbose=3)
if (!gsg$allOK) {
  Raw <- Raw[gsg$goodSample, gsg$goodGenes]
}
nrow(Raw)
ncol(Raw)

library(flashClust)
#datFT <- readRDS("Removed_Outlier_cut3_Volume_Change.rds")
powers <- c(c(1:10), seq(from=12, to=20, by=2))
nGene <- 90000
sft <- pickSoftThreshold(Raw, powerVector=powers, blockSize = nGene, verbose=5)

saveRDS(sft,file="./WGCNA/Power_information_t14_full_regress.rds")
########################################


sft <- readRDS('./WGCNA/Power_information_t14_full_regress.rds')
pdf(file="./WGCNA/Soft-threshold_t14_full_regress.pdf", width=9, height=5)
par(mfrow=c(1,2))   #Two plots
cex1=0.8 #具体看一下图

plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],xlab="Soft Threshold (power)", ylab="Scale Free Topology Model Fit, signed R^2", type="n", main=paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], labels=powers, cex=cex1, col="red")
abline(h=0.8, col="red")


powers <- c(c(1:10), seq(from=12, to=20, by=2))
nGene <- 90000
plot(sft$fitIndices[,1], sft$fitIndices[,5], xlab="Soft Threshold (power)", ylab= "Mean Connectivity", type="n", main= paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1, col="red")

dev.off()

################################


