
################
rm(list=ls())
setwd('/Users/dichen/Documents/Methylation_Running/List_ME/')
library(flashClust)
library(R.matlab)
r_substract <- read.csv('r_substract.csv')
head(r_substract)

library(cluster)
library(factoextra)
library(ggplot2)

##################


sampleTree <- flashClust(dist(r_substract,method = 'euclidean'), method="average")
saveRDS(sampleTree, file="./sampleTree.rds", compress=T)


plot(sampleTree, main= "WVCNA clusters clustering",labels = FALSE, hang = 1, xlab = 'hclust')
rect.hclust(sampleTree,k=18,border = 2:4)

Clust <- cutree(sampleTree,k=18) # 10% of 175 MEs is 18 clus

list(Clust)
tabulate(Clust)
length(unique(Clust))


pdf(file="Plot_Hierarchical_clustering_of_No_Merge_Network_Change.pdf", width=12, height=9)      
par(cex=0.4)
par(mar=c(0,4,2,0))
plot(sampleTree, main= "WVCNA clusters clustering",labels = FALSE, hang = 1, xlab = 'hclust')
rect.hclust(sampleTree,k=18,border = 2:4)
dev.off()

## save matrix_Module_Cluster
matrix_Module_Cluster<-cbind(c(1:175), Clust)
writeMat("../(5)matrix_Module_Cluster_percent_10.mat", matrix_Module_Cluster = matrix_Module_Cluster)
#####################################################################
####################### Heatmap
setwd('/Users/dichen/Documents/Methylation_Running/List_ME/')
library(pheatmap)

order_consensus_14<-as.double(sampleTree$order) #num of ME
length(order_consensus_14)


k<-order_consensus_14[2]
Clust[k] #num of class
tabulate(Clust)

## map to cluster num
num_cluster<-matrix(data=NA, nrow =175,ncol=1, byrow = FALSE, dimnames = NULL)
for (i in 1:175) {
  i
  k<-order_consensus_14[i]
  num_cluster[i]<-Clust[k] #num of class
}
unique(num_cluster)

## age 14
ME_consensus_14 <- read.csv('ME_consensus_14.csv')

ME_consensus_14_order=array(0,dim=dim(ME_consensus_14)) 
for (i in 1:length(order_consensus_14)){
  cmd<-order_consensus_14[i]
  ME_consensus_14_order[,i]<-ME_consensus_14[,cmd]
}

r_14_order<- cor(ME_consensus_14_order, method = "pearson", use = "complete.obs") 

#pheatmap(r_14_order, cluster_rows = T,cluster_cols = T,clustering_method = 'average')
pheatmap(r_14_order, cluster_rows = F,cluster_cols = F)

## age 19
ME_consensus_19 <- read.csv('ME_consensus_19.csv')

ME_consensus_19_order=array(0,dim=dim(ME_consensus_19)) 
for (i in 1:length(order_consensus_14)){
  cmd<-order_consensus_14[i]
  ME_consensus_19_order[,i]<-ME_consensus_19[,cmd]
}

#View(ME_consensus_19_order)
r_19_order<- cor(ME_consensus_19_order, method = "pearson", use = "complete.obs") 
pheatmap(r_19_order, cluster_rows = F,cluster_cols = F)

## Change
ME_consensus_substract_19_14 <- read.csv('ME_consensus_substract_19_14.csv')

ME_consensus_substract_19_14_order=array(0,dim=dim(ME_consensus_substract_19_14)) 
for (i in 1:length(order_consensus_14)){
  cmd<-order_consensus_14[i]
  ME_consensus_substract_19_14_order[,i]<-ME_consensus_substract_19_14[,cmd]
}

r_substract_19_14_order<- cor(ME_consensus_substract_19_14_order, method = "pearson", use = "complete.obs") 

pheatmap(r_substract_19_14_order, cluster_rows = T,cluster_cols = T,clustering_method = 'average')
pheatmap(r_substract_19_14_order, cluster_rows = F,cluster_cols = F)


## finally big figure # size=27*27

pheatmap(r_substract_19_14_order, cluster_rows = T,cutree_rows = 18,clustering_method = 'average')

rnames<-paste0("ME-", order_consensus_14)
cnames<-paste0("Gene", LETTERS[1:175])
r_data<-matrix(r_substract_19_14_order, nrow = 175,ncol = 175,dimnames=list(rnames, cnames))

pheatmap(r_data, cluster_rows = T,cutree_rows = 18,color = colorRampPalette(colors = c("#5EAD83","white","#9F558E"))(100),
         show_rownames = F,show_colnames = FALSE,cluster_cols = F,clustering_method = 'average')




