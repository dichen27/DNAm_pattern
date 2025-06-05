rm(list=ls())
setwd('/Users/dichen/Documents/Methylation_Running/')


#install.packages("R.matlab")
library(R.matlab)
num <- 372582

## t_14
t14_source<-readMat('data_regressed/t14_full_regress.mat')
t14<-t14_source$t14.full.regress
t14[1:3,1:3]
t_14<-t14[,-1]
t_14[1:3,1:3]
nrow(t_14)
ncol(t_14)

## t_19
t19_source<-readMat('data_regressed/t19_full_regress.mat')
t19<-t19_source$t19.full.regress
t19[1:3,1:3]
t_19<-t19[,-1]
t_19[1:3,1:3]
nrow(t_19)
ncol(t_19)

################## save rds
nSets <- 2 ## The number of dataset
multiExpr <- vector(mode = "list", length = nSets)
multiExpr[[1]] <- list(data = t_14)
multiExpr[[2]] <- list(data = t_19)
names(multiExpr) <- c("t14", "t19")

#saveRDS(multiExpr, file="data_regressed/Multi_Matrix_regress.rds",compress=T)

#################################
setwd('/Users/dichen/Documents/Methylation_Running/')

library(WGCNA)
allowWGCNAThreads()
options(stringsAsFactors=FALSE)

#Raw <- readRDS("data_regressed/Multi_Matrix_regress.rds")
Raw <- multiExpr
exprSize = checkSets(multiExpr)

nGene <- 40000 ## The number of voxels
bwnet <- blockwiseConsensusModules(
  Raw,maxBlockSize = nGene, corType = "pearson", power = 5,
  networkType = "unsigned", TOMType = "unsigned", TOMDenom = "min",
  saveIndividualTOMs = F, saveConsensusTOMs = T, deepSplit = 4,
  detectCutHeight = 0.995, minModuleSize = 50, pamRespectsDendro = FALSE,
  reassignThresholdPS = 0, mergeCutHeight = 0, numericLabels = T, 
  nThreads = 20, verbose = 3)
saveRDS(bwnet, file="Consensus_19_14_Network_Information_Power5.rds", compress=T)
######################

####################
## The option is the same as the step3_WVCNA_main

setwd('/Users/dichen/Documents/Methylation_Running/')

bwnet <- readRDS('Consensus_19_14_Network_Information_Power5.rds')
Origin <- as.data.frame(bwnet$colors)
write.table(Origin, file="List_Of_Methylation_Change_Power5_Min50_Merge0_UnsignedTOM_Consensus_19_14.txt", quote=F, sep="\t", col.names=F, row.names=F
###########