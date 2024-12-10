rm(list=ls())
library(R.matlab) 
library(xlsx)
library(pheatmap)
setwd("/Users/dichen/Documents/ADNI_Running/")

##########################
## load the ADNI data-606
load("./data_ADNI2_606/data_ADNI_606.rda")

## cov
load("./data_ADNI2_606/cov_ADNI_methy_606.rda")
cov_ADNI_methy_606[1:3,]
write.csv(cov_ADNI_methy_606, file = "./data_ADNI2_606/cov_ADNI_methy_606.csv") # save for matlab

## ID
load("./data_ADNI2_606/ID_long_ADNI_606.rda")
load("./data_ADNI2_606/ID_short_ADNI_606.rda")
load("./data_ADNI2_606/cpg_name_ADNI.rda")

## grouplabel
load("./data_ADNI2_606/grouplabel_606.rda")
grouplabel_606$grouplabel
write.csv(grouplabel_606, file = "./data_ADNI2_606/grouplabel_606.csv") # save for matlab
##################### load IMAGEN result
bwnet <- readRDS('./Cluster18_IMAGEN/Consensus_19_14_Network_Information_Power5.rds')
Origin <- as.data.frame(bwnet$colors)
index_CpG<-bwnet$colors


load('./Cluster18_IMAGEN/t14_506_ori.rda')
t14_2[1:3,1:3]
all_cpg_IMAGEN<-colnames(t14_2)
all_cpg_IMAGEN<-all_cpg_IMAGEN[2:length(all_cpg_IMAGEN)] # dele the "ID"

################### mapping to ME-175
ADNI_ME175<-matrix(data=NA, nrow = length(ID_short_ADNI_2_unique), ncol =175, byrow = FALSE, dimnames = NULL)
for (k in 1:175){
 
  print(k)
  
  ID_CpGs_k<-all_cpg_IMAGEN[index_CpG==k]
  
  data_k<-matrix(data=NA, nrow = length(ID_short_ADNI_2_unique), ncol =length(ID_CpGs_k), byrow = FALSE, dimnames = NULL)
  for (i in 1:length(ID_CpGs_k)){
    
    lin<-data_ADNI_606[,cpg_name_ADNI==ID_CpGs_k[i]]
    data_k[,i]<-lin
  }
  data_k_mean<-apply(data_k,1,mean)  # mean of the each row
  ADNI_ME175[,k]<-data_k_mean
}

ADNI_ME175<-cbind(ID_long_ADNI_2_unique,ID_short_ADNI_2_unique,ADNI_ME175)
write.csv(ADNI_ME175, file = "ADNI_ME175.csv")

ADNI_ME175[1:3,1:3]
dim(ADNI_ME175)
################### mapping to Cluster-18

ADNI_cluster18<-matrix(data=NA, nrow = length(ID_short_ADNI_2_unique), ncol =18, byrow = FALSE, dimnames = NULL)
for (k in 1:18){
  print(k)
  cluster_dir<-paste("./Cluster18_IMAGEN/CLuster_CpGs_k",k,".rds",sep = "")
  ID_CpGs_k <- readRDS(cluster_dir)
  data_k<-matrix(data=NA, nrow = length(ID_short_ADNI_2_unique), ncol =length(ID_CpGs_k), byrow = FALSE, dimnames = NULL)
  for (i in 1:length(ID_CpGs_k)){
    
    lin<-data_ADNI_606[,cpg_name_ADNI==ID_CpGs_k[i]]
    data_k[,i]<-lin
  }
  data_k_mean<-apply(data_k,1,mean)  # mean of the each row
  ADNI_cluster18[,k]<-data_k_mean
}

ADNI_cluster18_ID<-cbind(ID_short_ADNI_2_unique,ADNI_cluster18)
write.csv(ADNI_cluster18_ID, file = "ADNI_cluster18_ID.csv")

dim(ADNI_cluster18_ID)
ADNI_cluster18_ID[1:3,1:3]
########################################################## 
## Mantel test

## load IMAGEN final pattern
mat_data <- readMat("../Methylation_Running/(Reply1-8)pattern_longitudinal_cluster18_regressed.mat")
r_IMAGEN_change<- mat_data$pattern.longitudinal.cluster18.regressed



## Full AD data
r_ADNI_all_cluster18<- cor(ADNI_cluster18, method = "pearson", use = "complete.obs") 
pheatmap(r_ADNI_all_cluster18, cluster_rows = F,cluster_cols = F)
pheatmap(r_ADNI_all_cluster18, cluster_rows = F,cluster_cols = F,color = colorRampPalette(colors = c("#5EAD83","white","#9F558E"))(100))


library(vegan)
mantel(xdis = r_IMAGEN_change, ydis = r_ADNI_all_cluster18, method = "pearson", permutations = 10000, na.rm = TRUE)

## read grouping data
grouping<-readMat('./(3)ADNI_ID_grouping.mat')

## AD
ID_AD<-grouping$ID.AD
length(ID_AD)
ADNI_cluster18_AD<-matrix(data=NA, nrow =length(ID_AD), ncol =18, byrow = FALSE, dimnames = NULL)
for (i in 1:length(ID_AD)){
  ADNI_cluster18_AD[i,]<-ADNI_cluster18_ID[ADNI_cluster18_ID[,1]==ID_AD[i],-1]
}
r_ADNI_cluster18_AD<- cor(ADNI_cluster18_AD, method = "pearson", use = "complete.obs") 

mantel(xdis = r_IMAGEN_change, ydis = r_ADNI_cluster18_AD, method = "pearson", permutations = 10000, na.rm = TRUE)

## MCI
ID_MCI<-grouping$ID.MCI
length(ID_MCI)
ADNI_cluster18_MCI<-matrix(data=NA, nrow =length(ID_MCI), ncol =18, byrow = FALSE, dimnames = NULL)
for (i in 1:length(ID_MCI)){
  ADNI_cluster18_MCI[i,]<-ADNI_cluster18_ID[ADNI_cluster18_ID[,1]==ID_MCI[i],-1]
}
r_ADNI_cluster18_MCI<- cor(ADNI_cluster18_MCI, method = "pearson", use = "complete.obs") 

mantel(xdis = r_IMAGEN_change, ydis = r_ADNI_cluster18_MCI, method = "pearson", permutations = 10000, na.rm = TRUE)

## HC
ID_HC<-grouping$ID.HC
length(ID_HC)
ADNI_cluster18_HC<-matrix(data=NA, nrow =length(ID_HC), ncol =18, byrow = FALSE, dimnames = NULL)
for (i in 1:length(ID_HC)){
  ADNI_cluster18_HC[i,]<-ADNI_cluster18_ID[ADNI_cluster18_ID[,1]==ID_HC[i],-1]
}
r_ADNI_cluster18_HC<- cor(ADNI_cluster18_HC, method = "pearson", use = "complete.obs") 

mantel(xdis = r_IMAGEN_change, ydis = r_ADNI_cluster18_HC, method = "pearson", permutations = 10000, na.rm = TRUE)

########################################################## 
## plot the heatmap
## order the cluster
num_order<-c(4,1,13,12,6,11,18,8,5,2,14,3,15,10,9,7,16,17)

ADNI_cluster18_order<-matrix(data=NA, nrow =606, ncol =18, byrow = FALSE, dimnames = NULL)
for (i in 1:18){
  k<-num_order[i]
  ADNI_cluster18_order[,i]<-ADNI_cluster18[,k]
}
ADNI_cluster18_order_ID<-cbind(ID_short_ADNI_2_unique,ADNI_cluster18_order)

## Full PPMI data
r_ADNI_all_cluster18_order<- cor(ADNI_cluster18_order, method = "pearson", use = "complete.obs") 
pheatmap(r_ADNI_all_cluster18_order, cluster_rows = F,cluster_cols = F)

rnames<-paste0("Cluster-", num_order)
cnames<-paste0("Cluster-", num_order)
r_data<-matrix(r_ADNI_all_cluster18_order, nrow = 18,ncol = 18,dimnames=list(rnames, cnames))
pheatmap(r_data, cluster_rows = F,cluster_cols = F,show_rownames = T,show_colnames = F,border_color = "white",color = colorRampPalette(colors = c("#5EAD83","white","#9F558E"))(100))

## AD
ID_AD<-grouping$ID.AD
length(ID_AD)
ADNI_cluster18_order_AD<-matrix(data=NA, nrow =length(ID_AD), ncol =18, byrow = FALSE, dimnames = NULL)
for (i in 1:length(ID_AD)){
  ADNI_cluster18_order_AD[i,]<-ADNI_cluster18_order_ID[ADNI_cluster18_order_ID[,1]==ID_AD[i],-1]
}
r_ADNI_cluster18_order_AD<- cor(ADNI_cluster18_order_AD, method = "pearson", use = "complete.obs") 

r_data<-matrix(r_ADNI_cluster18_order_AD, nrow = 18,ncol = 18,dimnames=list(rnames, cnames))
pheatmap(r_data, cluster_rows = F,cluster_cols = F,show_rownames = T,show_colnames = F,border_color = "white",color = colorRampPalette(colors = c("#5EAD83","white","#9F558E"))(100))

## MCI
ID_MCI<-grouping$ID.MCI
length(ID_MCI)

ADNI_cluster18_order_MCI<-matrix(data=NA, nrow =length(ID_MCI), ncol =18, byrow = FALSE, dimnames = NULL)
for (i in 1:length(ID_MCI)){
  ADNI_cluster18_order_MCI[i,]<-ADNI_cluster18_order_ID[ADNI_cluster18_order_ID[,1]==ID_MCI[i],-1]
}
r_ADNI_cluster18_order_MCI<- cor(ADNI_cluster18_order_MCI, method = "pearson", use = "complete.obs") 

r_data<-matrix(r_ADNI_cluster18_order_MCI, nrow = 18,ncol = 18,dimnames=list(rnames, cnames))
pheatmap(r_data, cluster_rows = F,cluster_cols = F,show_rownames = T,show_colnames = F,border_color = "white",color = colorRampPalette(colors = c("#5EAD83","white","#9F558E"))(100))

## HC
ID_HC<-grouping$ID.HC
length(ID_HC)

ADNI_cluster18_order_HC<-matrix(data=NA, nrow =length(ID_HC), ncol =18, byrow = FALSE, dimnames = NULL)
for (i in 1:length(ID_HC)){
  ADNI_cluster18_order_HC[i,]<-ADNI_cluster18_order_ID[ADNI_cluster18_order_ID[,1]==ID_HC[i],-1]
}
r_ADNI_cluster18_order_HC<- cor(ADNI_cluster18_order_HC, method = "pearson", use = "complete.obs") 

r_data<-matrix(r_ADNI_cluster18_order_HC, nrow = 18,ncol = 18,dimnames=list(rnames, cnames))
pheatmap(r_data, cluster_rows = F,cluster_cols = F,show_rownames = T,show_colnames = F,border_color = "white",color = colorRampPalette(colors = c("#5EAD83","white","#9F558E"))(100))


##################### finished









