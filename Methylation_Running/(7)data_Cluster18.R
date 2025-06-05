rm(list=ls())
setwd('/Users/dichen/Documents/Methylation_Running/')

library(missMethyl)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(limma)


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

################### age-14

t14_matrix<-t14_2[,2:ncol(t14_2)]
dim(t14_matrix)
t14_matrix[1:3,1:3]

data_14_origin_cluster18<-matrix(data=NA, nrow = 506, ncol =18, byrow = FALSE, dimnames = NULL)
for (k in 1:18){
  print(k)
cluster_dir<-paste("./Cluster18/CLuster_CpGs_k",k,".rds",sep = "")
ID_CpGs_k <- readRDS(cluster_dir)

data_k<-matrix(data=NA, nrow = 506, ncol =length(ID_CpGs_k), byrow = FALSE, dimnames = NULL)
for (i in 1:length(ID_CpGs_k)){
lin<-t14_matrix[,colnames(t14_matrix)==ID_CpGs_k[i]]
data_k[,i]<-lin
}

data_k_mean<-apply(data_k,1,mean)  # mean of the each row
data_14_origin_cluster18[,k]<-data_k_mean
}

write.csv(data_14_origin_cluster18, file = "data_14_origin_cluster18.csv")

###################### age-19
load('./IMAGEN_methylation_data_origin/t19_506_ori.rda')
t19_2[1:3,1:3]
dim(t19_2)
t19_matrix<-t19_2[,2:ncol(t19_2)]
dim(t19_matrix)
t19_matrix[1:3,1:3]


data_19_origin_cluster18<-matrix(data=NA, nrow = 506, ncol =18, byrow = FALSE, dimnames = NULL)
for (k in 1:18){
  print(k)
  cluster_dir<-paste("./Cluster18/CLuster_CpGs_k",k,".rds",sep = "")
  ID_CpGs_k <- readRDS(cluster_dir)
  
  data_k<-matrix(data=NA, nrow = 506, ncol =length(ID_CpGs_k), byrow = FALSE, dimnames = NULL)
  for (i in 1:length(ID_CpGs_k)){
    lin<-t19_matrix[,colnames(t19_matrix)==ID_CpGs_k[i]]
    data_k[,i]<-lin
  }
  
  data_k_mean<-apply(data_k,1,mean)  # mean of the each row
  data_19_origin_cluster18[,k]<-data_k_mean
}

write.csv(data_19_origin_cluster18, file = "data_19_origin_cluster18.csv")
################### data change
data_change_origin_cluster18<-data_19_origin_cluster18-data_14_origin_cluster18
write.csv(data_change_origin_cluster18, file = "data_change_origin_cluster18.csv")










