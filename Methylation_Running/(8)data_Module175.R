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

###################### load the index of 175 MEs
bwnet <- readRDS('./Consensus_19_14_Network_Information_Power5.rds')
Origin <- as.data.frame(bwnet$colors)

# save for Li-Ao
saveRDS(Origin, file = "./Module176_index.rds")
# kkk <- readRDS('./Module176_index.rds')
index_CpG<-bwnet$colors


length(unique(index_CpG))
sum(index_CpG==3)
sig_cpg_cmd<-all_cpg_IMAGEN[index_CpG==3]
length(sig_cpg_cmd)
###################
##  save Module175 result
for (k in 1:175){

print(k)
ID_CpGs_k<-all_cpg_IMAGEN[index_CpG==k]
ID_CpGs_k[1:6]

Module_dir<-paste("./Module175/Module_CpGs_k",k,".rds",sep = "")
saveRDS(ID_CpGs_k, Module_dir)
}


################### age-14
t14_matrix<-t14_2[,2:ncol(t14_2)]
dim(t14_matrix)
t14_matrix[1:3,1:3]

data_14_origin_ME175<-matrix(data=NA, nrow = 506, ncol =175, byrow = FALSE, dimnames = NULL)
for (k in 1:175){
  print(k)

ID_CpGs_k<-all_cpg_IMAGEN[index_CpG==k]
data_k<-matrix(data=NA, nrow = 506, ncol =length(ID_CpGs_k), byrow = FALSE, dimnames = NULL)
for (i in 1:length(ID_CpGs_k)){
lin<-t14_matrix[,colnames(t14_matrix)==ID_CpGs_k[i]]
data_k[,i]<-lin
}

data_k_mean<-apply(data_k,1,mean)  # mean of the each row
data_14_origin_ME175[,k]<-data_k_mean
}

write.csv(data_14_origin_ME175, file = "data_14_origin_ME175.csv")

###################### age-19
load('./IMAGEN_methylation_data_origin/t19_506_ori.rda')
t19_2[1:3,1:3]
dim(t19_2)
t19_matrix<-t19_2[,2:ncol(t19_2)]
dim(t19_matrix)
t19_matrix[1:3,1:3]


data_19_origin_ME175<-matrix(data=NA, nrow = 506, ncol =175, byrow = FALSE, dimnames = NULL)
for (k in 1:175){
  print(k)
  
  ID_CpGs_k<-all_cpg_IMAGEN[index_CpG==k]
  data_k<-matrix(data=NA, nrow = 506, ncol =length(ID_CpGs_k), byrow = FALSE, dimnames = NULL)
  for (i in 1:length(ID_CpGs_k)){
    lin<-t19_matrix[,colnames(t19_matrix)==ID_CpGs_k[i]]
    data_k[,i]<-lin
  }
  
  data_k_mean<-apply(data_k,1,mean)  # mean of the each row
  data_19_origin_ME175[,k]<-data_k_mean
}

write.csv(data_19_origin_ME175, file = "data_19_origin_ME175.csv")
################ 













