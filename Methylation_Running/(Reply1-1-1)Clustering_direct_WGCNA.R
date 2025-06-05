
rm(list=ls())
setwd('/Users/dichen/Documents/Methylation_Running/')
#setwd('F:/Methylation_Running/')

library(R.matlab)

num <- 372582
sub <- 506
mat_data <- readMat("./List_ME/substract_19_14_all_CpGs.mat")
substract_19_14_all_CpGs<-mat_data$substract.19.14.all.CpGs

Raw<-substract_19_14_all_CpGs[,-1]




library(WGCNA)
allowWGCNAThreads()
options(stringsAsFactors=FALSE)
gsg <- goodSamplesGenes(as.matrix(Raw), verbose=3)
if (!gsg$allOK) {
  Raw <- Raw[gsg$goodSample, gsg$goodGenes]
}



library(flashClust)

#datFT <- readRDS("Removed_Outlier_cut3_Volume_Change.rds")
powers <- c(c(1:10), seq(from=12, to=20, by=2))
nGene <- 90000
sft <- pickSoftThreshold(Raw, powerVector=powers, blockSize = nGene, verbose=5)

saveRDS(sft,file="./WGCNA/Power_information_substract_19_14_all_CpGs.rds")
########################################



##############################

sft <- readRDS('./WGCNA/Power_information_substract_19_14_all_CpGs.rds')
pdf(file="./WGCNA/Soft-threshold_substract_19_14_all_CpGs.pdf", width=9, height=5)
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

#######################

library(WGCNA)
allowWGCNAThreads()
options(stringsAsFactors=FALSE)
#datFT <- readRDS("Removed_Outlier_cut3_Volume_Change.rds")
nGene <- 40000


bwnet <- blockwiseModules(Raw, networkType="unsigned", TOMType="unsigned", power=3, 
                          maxBlockSize = nGene, detectCutHeight=0.995, minModuleSize=50, 
                          reassignThreshold=0, deepSplit=4, mergeCutHeight=0, numericLabels=T, 
                          saveTOMs=T, nThreads=20, verbose=3) 
saveRDS(bwnet, file="./WGCNA/substract_19_14_all_CpGs_Power3.rds", compress=T)
#######################


#######################
bwnet <- readRDS('./WGCNA/substract_19_14_all_CpGs_Power3.rds')
Origin <- as.data.frame(bwnet$colors)

Origin$`bwnet$colors`
length(tabulate(Origin$`bwnet$colors`))

write.table(Origin, file="./WGCNA/List_Of_Methylation_Change_Power3_Min50_Merge0_UnsignedTOM_substract_19_14_all_CpGs.txt", quote=F, sep="\t", col.names=F, row.names=F)

ME <- as.data.frame(bwnet$MEs)
write.table(ME, file="./WGCNA/ME_Methylation_Change_Power3_Min50_Merge0_UnsignedTOM_substract_19_14_all_CpGs.txt", quote=F, sep="\t", col.names=T, row.names=F)
##############################
##########################
library(flashClust)

r <- cor(ME, method = "pearson", use = "complete.obs")         

sampleTree <- flashClust(dist(r), method="average")

pdf(file="./WGCNA/Plot_Hierarchical_clustering_of_Power3_No_Merge_Network_substract_19_14_all_CpGs.pdf", width=12, height=9)      
par(cex=0.4)
par(mar=c(0,4,2,0))
plot(sampleTree, main= "WVCNA clusters clustering", sub="", xlab="", cex.lab=0.5, cex.axis=0.5, cex.main=2)
abline (h=1000, col="red")
dev.off()

library(pheatmap)
pheatmap(r, cluster_rows = T,cluster_cols = T,clustering_method = 'average')
