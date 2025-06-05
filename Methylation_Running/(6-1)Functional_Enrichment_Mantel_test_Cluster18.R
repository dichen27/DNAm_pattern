rm(list=ls())
setwd('/Users/dichen/Documents/Methylation_Running/')

library(missMethyl)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(limma)
library(ggplot2)
library(pheatmap)
library(vegan)
###################### missMehtyl example data
ann <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)

# All CpG sites tested
allcpgs <- rownames(ann)
allcpgs[1]

##################### load all IMAGEN CpGs
# age 14
load('./IMAGEN_methylation_data_origin/t14_506_ori.rda')
t14_2[1:3,1:3]


all_cpg_IMAGEN<-colnames(t14_2)
#all_cpg_IMAGEN<-colnames(all_cpg_IMAGEN)

all_cpg_IMAGEN[1:3]
all_cpg_IMAGEN<-all_cpg_IMAGEN[2:length(all_cpg_IMAGEN)] # dele the "ID"

length(all_cpg_IMAGEN)
all_cpg_IMAGEN[1:3]

################### load the consensus result
bwnet <- readRDS('Consensus_19_14_Network_Information_Power5.rds')
Origin <- as.data.frame(bwnet$colors)
index_CpG<-bwnet$colors


ME_consensus_14 <-read.csv(paste(getwd(),"/List_ME/ME_consensus_14.csv",sep = ""))
ME_consensus_14<-data.matrix(ME_consensus_14)
dim(ME_consensus_14)

################### load second cluster result

library(flashClust)

r_substract <- read.csv(paste(getwd(),"/List_ME/r_substract.csv",sep = ""))
r_substract[1:3,1:3]

sampleTree <- flashClust(dist(r_substract,method = 'euclidean'), method="average")

plot(sampleTree, main= "WVCNA clusters clustering", sub="", xlab="", cex.lab=0.5, cex.axis=0.5, cex.main=2)
#plot(sampleTree, main= "WVCNA clusters clustering",labels = FALSE, hang = 1, xlab = 'hclust')
rect.hclust(sampleTree,k=18,border = 2:4)


Clust <- cutree(sampleTree,k=18) # 10% of 175 MEs is 18 clus
write.csv(Clust, file = "Cluster_index_IMAGEN.csv")

order_consensus_14<-as.double(sampleTree$order)

tabulate(Clust)

Num_Cluster_18<-tabulate(Clust)

order_consensus_14
num_175<-c(1:175)

############ align the data to 506*18


## load change pattern
IMAGEN_change_cluster18<-read.csv('./data_change_origin_cluster18.csv')
IMAGEN_change_cluster18<-IMAGEN_change_cluster18[,-1]
r_IMAGEN_change<- cor(IMAGEN_change_cluster18, method = "pearson", use = "complete.obs") 
pheatmap(r_IMAGEN_change, cluster_rows = F,cluster_cols = F)


# load the new data_14
data_14_cluster18<-read.csv('./data_14_origin_cluster18.csv') 
data_14_cluster18<-data_14_cluster18[,-1]

r_14_cluster18<- cor(data_14_cluster18, method = "pearson", use = "complete.obs") 
pheatmap(r_14_cluster18, cluster_rows = F,cluster_cols = F)

## mantel test
kkk<-mantel(xdis = r_IMAGEN_change, ydis = r_14_cluster18, method = "pearson", permutations = 10000, na.rm = TRUE)
kkk$statistic
## order by change pattern
num_order<-c(4,1,13,12,6,11,18,8,5,2,14,3,15,10,9,7,16,17)

data_14_cluster18_order<-matrix(data=NA, nrow =506, ncol =18, byrow = FALSE, dimnames = NULL)
for (i in 1:18){
  k<-num_order[i]
  data_14_cluster18_order[,i]<-data_14_cluster18[,k]
}

## heatmap
r_IMAGEN_14_order<- cor(data_14_cluster18_order, method = "pearson", use = "complete.obs") 
pheatmap(r_IMAGEN_14_order, cluster_rows = F,cluster_cols = F)

rnames<-paste0("Cluster-", num_order)
cnames<-paste0("Cluster-", num_order)
r_data<-matrix(r_IMAGEN_14_order, nrow = 18,ncol = 18,dimnames=list(rnames, cnames))
pheatmap(r_data, cluster_rows = F,cluster_cols = F,show_rownames = T,show_colnames = F,border_color = "white",color = colorRampPalette(colors = c("#5EAD83","white","#9F558E"))(100))



########################################
## age 19
# load the new data_19
data_19_cluster18<-read.csv('./data_19_origin_cluster18.csv') 
data_19_cluster18<-data_19_cluster18[,-1]


## mantel test
r_19_cluster18<- cor(data_19_cluster18, method = "pearson", use = "complete.obs") 
pheatmap(r_19_cluster18, cluster_rows = F,cluster_cols = F)


mantel(xdis = r_IMAGEN_change, ydis = r_19_cluster18, method = "pearson", permutations = 10000, na.rm = TRUE)

## order by change pattern
num_order<-c(4,1,13,12,6,11,18,8,5,2,14,3,15,10,9,7,16,17)

data_19_cluster18_order<-matrix(data=NA, nrow =506, ncol =18, byrow = FALSE, dimnames = NULL)
for (i in 1:18){
  k<-num_order[i]
  data_19_cluster18_order[,i]<-data_19_cluster18[,k]
}


## heatmap
r_IMAGEN_19_order<- cor(data_19_cluster18_order, method = "pearson", use = "complete.obs") 
pheatmap(r_IMAGEN_19_order, cluster_rows = F,cluster_cols = F)

rnames<-paste0("Cluster-", num_order)
cnames<-paste0("Cluster-", num_order)
r_data<-matrix(r_IMAGEN_19_order, nrow = 18,ncol = 18,dimnames=list(rnames, cnames))
pheatmap(r_data, cluster_rows = F,cluster_cols = F,show_rownames = T,show_colnames = F,border_color = "white",color = colorRampPalette(colors = c("#5EAD83","white","#9F558E"))(100))


#################### Enrichment
num<-matrix(1:175,175,1)
num_CpG<-matrix(data=NA, nrow = length(unique(Clust)), ncol =1, byrow = FALSE, dimnames = NULL)

for (k in 1:length(unique(Clust))) {
#k<-8

print(paste("k=",k,sep = "")) # important function paste
cluster_k <- num[Clust==k]
length(cluster_k)

num_CpG_k<-0
sig_cpg_IMAGEN<-""
for (i in 1:length(cluster_k)) {

cmd<-cluster_k[i]


sig_cpg_cmd<-all_cpg_IMAGEN[index_CpG==cmd]
sig_cpg_IMAGEN<-c(sig_cpg_IMAGEN,sig_cpg_cmd)

}

sig_cpg_IMAGEN<-sig_cpg_IMAGEN[-1] # dele the ""



num_CpG[k]<-length(sig_cpg_IMAGEN)
print(paste("num of CpGs:",length(sig_cpg_IMAGEN)))
sig_cpg_IMAGEN[1]

### save the CpGs of Cluster
cluster5_dir<-paste(getwd(),"/Cluster18/CLuster_CpGs_k",k,".rds",sep = "")
saveRDS(sig_cpg_IMAGEN, file=cluster5_dir,compress=T)

sig_cpg_IMAGEN_k<- readRDS(cluster5_dir)
#length(sig_cpg_IMAGEN_k)

########################
# GO testing with prior probabilities taken into account
# Plot of bias due to differing numbers of CpG sites per gene
gst <- gometh(sig.cpg = sig_cpg_IMAGEN, all.cpg = all_cpg_IMAGEN, collection = "GO", 
              plot.bias = TRUE, prior.prob = TRUE, anno = ann)
# Total number of GO categories significant at 5% FDR
table(gst$FDR<0.05)
# Table of top GO results
topGSA(gst,n=10)

#### save GO
gst_dir<-paste(getwd(),"/Enrichment_Second_Round_Cluster18/Enrichment_Second_Round_GO_k",k,".rds",sep = "")
saveRDS(gst, file=gst_dir,compress=T)


# KEGG testing
kegg <- gometh(sig.cpg = sig_cpg_IMAGEN, all.cpg = all_cpg_IMAGEN, collection = "KEGG", 
               prior.prob=TRUE, anno = ann)
# FDR
table(kegg$FDR<0.05)
# Table of top KEGG results
topGSA(kegg,n=10)


#### save KEGG
kegg_dir<-paste(getwd(),"/Enrichment_Second_Round_Cluster18/Enrichment_Second_Round_KEGG_k",k,".rds",sep = "")
saveRDS(kegg, file=kegg_dir,compress=T)
}

########################







