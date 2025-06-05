
rm(list=ls())
library(R.matlab) 
library(pheatmap)
setwd("/Users/dichen/Documents/PPMI_Running/")


data_dir="/Users/dichen/Documents/PPMI_Running/PPMI/preprocess/"
##########################
## load the PPMI data

load(paste(data_dir,"Beta_QCed.rda",sep=""))
data_PPMI<-t(beta)
dim(data_PPMI)

## ID
load(paste(data_dir,"ID_short_final_QCed.rda",sep=""))
length(ID_short_final)

## cellcount
load(paste(data_dir,"cellcount_QCed.rda",sep=""))
dim(cellcount)

## PC-4
load(paste(data_dir,"fast_svd_QCed.rda",sep=""))
PCs<-ss$v[,1:4]
dim(PCs)

##   bind the cov
cov_PPMI_methy<-cbind(ID_short_final,cellcount,PCs)
write.csv(cov_PPMI_methy, file = "cov_PPMI_methy.csv")

################### mapping to cluster
PPMI_cluster18<-matrix(data=NA, nrow = length(ID_short_final), ncol =18, byrow = FALSE, dimnames = NULL)
for (k in 1:18){
  print(k)
  cluster_dir<-paste("./Cluster18_IMAGEN/CLuster_CpGs_k",k,".rds",sep = "")
  ID_CpGs_k <- readRDS(cluster_dir)
  data_k<-matrix(data=NA, nrow = length(ID_short_final), ncol =length(ID_CpGs_k), byrow = FALSE, dimnames = NULL)
  for (i in 1:length(ID_CpGs_k)){
    
    lin<-data_PPMI[,colnames(data_PPMI)==ID_CpGs_k[i]]
    data_k[,i]<-lin
  }
  data_k_mean<-apply(data_k,1,mean)  # mean of the each row
  PPMI_cluster18[,k]<-data_k_mean
}

PPMI_cluster18<-cbind(ID_short_final,PPMI_cluster18)
write.csv(PPMI_cluster18, file = "PPMI_cluster18.csv")

dim(PPMI_cluster18)
PPMI_cluster18[1:6,1:6]
########################################################## 
## Mantel test
## load IMAGEN regress pattern
IMAGEN_change_cluster18<-read.csv('../Methylation_Running/data_change_origin_cluster18.csv')
mat_data <- readMat("../Methylation_Running/(Reply1-8)pattern_longitudinal_cluster18_regressed.mat")
r_IMAGEN_change<- mat_data$pattern.longitudinal.cluster18.regressed


## Full PPMI data
PPMI_all_cluster18<-PPMI_cluster18[,-1]
r_PPMI_all_cluster18<- cor(PPMI_all_cluster18, method = "pearson", use = "complete.obs") 
pheatmap(r_PPMI_all_cluster18, cluster_rows = F,cluster_cols = F)

library(vegan)
mantel(xdis = r_IMAGEN_change, ydis = r_PPMI_all_cluster18, method = "pearson", permutations = 10000, na.rm = TRUE)

## read grouping data
grouping<-readMat('./(2-2)ID_grouping.mat')


## PD
ID_AD<-grouping$ID.AD
length(ID_AD)
PPMI_cluster18_AD<-matrix(data=NA, nrow =length(ID_AD), ncol =18, byrow = FALSE, dimnames = NULL)
for (i in 1:length(ID_AD)){
PPMI_cluster18_AD[i,]<-PPMI_cluster18[PPMI_cluster18[,1]==ID_AD[i],-1]
}
r_PPMI_cluster18_AD<- cor(PPMI_cluster18_AD, method = "pearson", use = "complete.obs") 

mantel(xdis = r_IMAGEN_change, ydis = r_PPMI_cluster18_AD, method = "pearson", permutations = 10000, na.rm = TRUE)

## SWEDD
ID_SWEDD<-grouping$ID.SWEDD
length(ID_SWEDD)
PPMI_cluster18_SWEDD<-matrix(data=NA, nrow =length(ID_SWEDD), ncol =18, byrow = FALSE, dimnames = NULL)
for (i in 1:length(ID_SWEDD)){
  PPMI_cluster18_SWEDD[i,]<-PPMI_cluster18[PPMI_cluster18[,1]==ID_SWEDD[i],-1]
}
r_PPMI_cluster18_SWEDD<- cor(PPMI_cluster18_SWEDD, method = "pearson", use = "complete.obs") 
mantel(xdis = r_IMAGEN_change, ydis = r_PPMI_cluster18_SWEDD, method = "pearson", permutations = 10000, na.rm = TRUE)

## HC
ID_HC<-grouping$ID.HC
length(ID_HC)
PPMI_cluster18_HC<-matrix(data=NA, nrow =length(ID_HC), ncol =18, byrow = FALSE, dimnames = NULL)
for (i in 1:length(ID_HC)){
  PPMI_cluster18_HC[i,]<-PPMI_cluster18[PPMI_cluster18[,1]==ID_HC[i],-1]
}
r_PPMI_cluster18_HC<- cor(PPMI_cluster18_HC, method = "pearson", use = "complete.obs") 
mantel(xdis = r_IMAGEN_change, ydis = r_PPMI_cluster18_HC, method = "pearson", permutations = 10000, na.rm = TRUE)
########################################################## 
## plot the heatmap
## order the cluster
num_order<-c(4,1,13,12,6,11,18,8,5,2,14,3,15,10,9,7,16,17)

PPMI_cluster18_order<-matrix(data=NA, nrow =518, ncol =18, byrow = FALSE, dimnames = NULL)
for (i in 1:18){
  k<-num_order[i]
  PPMI_cluster18_order[,i]<-PPMI_all_cluster18[,k]
}

PPMI_cluster18_order_ID<-cbind(PPMI_cluster18[,1],PPMI_cluster18_order)

## Full PPMI data
r_PPMI_all_cluster18_order<- cor(PPMI_cluster18_order, method = "pearson", use = "complete.obs") 
pheatmap(r_PPMI_all_cluster18_order, cluster_rows = F,cluster_cols = F)

rnames<-paste0("Cluster-", num_order)
cnames<-paste0("Cluster-", num_order)
r_data<-matrix(r_PPMI_all_cluster18_order, nrow = 18,ncol = 18,dimnames=list(rnames, cnames))
pheatmap(r_data, cluster_rows = F,cluster_cols = F,show_rownames = T,show_colnames = F,border_color = "white",color = colorRampPalette(colors = c("#5EAD83","white","#9F558E"))(100))


## PD
ID_AD<-grouping$ID.AD
length(ID_AD)
PPMI_cluster18_order_AD<-matrix(data=NA, nrow =length(ID_AD), ncol =18, byrow = FALSE, dimnames = NULL)
for (i in 1:length(ID_AD)){
  PPMI_cluster18_order_AD[i,]<-PPMI_cluster18_order_ID[PPMI_cluster18_order_ID[,1]==ID_AD[i],-1]
}
r_PPMI_cluster18_order_AD<- cor(PPMI_cluster18_order_AD, method = "pearson", use = "complete.obs") 

r_data<-matrix(r_PPMI_cluster18_order_AD, nrow = 18,ncol = 18,dimnames=list(rnames, cnames))
pheatmap(r_data, cluster_rows = F,cluster_cols = F,show_rownames = T,show_colnames = F,border_color = "white",color = colorRampPalette(colors = c("#5EAD83","white","#9F558E"))(100))


## SWEDD
ID_SWEDD<-grouping$ID.SWEDD
length(ID_SWEDD)
PPMI_cluster18_order_SWEDD<-matrix(data=NA, nrow =length(ID_SWEDD), ncol =18, byrow = FALSE, dimnames = NULL)
for (i in 1:length(ID_SWEDD)){
  PPMI_cluster18_order_SWEDD[i,]<-PPMI_cluster18_order_ID[PPMI_cluster18_order_ID[,1]==ID_SWEDD[i],-1]
}
r_PPMI_cluster18_order_SWEDD<- cor(PPMI_cluster18_order_SWEDD, method = "pearson", use = "complete.obs") 

r_data<-matrix(r_PPMI_cluster18_order_SWEDD, nrow = 18,ncol = 18,dimnames=list(rnames, cnames))
pheatmap(r_data, cluster_rows = F,cluster_cols = F,show_rownames = T,show_colnames = F,border_color = "white",color = colorRampPalette(colors = c("#5EAD83","white","#9F558E"))(100))


## HC
ID_HC<-grouping$ID.HC
length(ID_HC)
PPMI_cluster18_order_HC<-matrix(data=NA, nrow =length(ID_HC), ncol =18, byrow = FALSE, dimnames = NULL)
for (i in 1:length(ID_HC)){
  PPMI_cluster18_order_HC[i,]<-PPMI_cluster18_order_ID[PPMI_cluster18_order_ID[,1]==ID_HC[i],-1]
}
r_PPMI_cluster18_order_HC<- cor(PPMI_cluster18_order_HC, method = "pearson", use = "complete.obs") 

r_data<-matrix(r_PPMI_cluster18_order_HC, nrow = 18,ncol = 18,dimnames=list(rnames, cnames))
pheatmap(r_data, cluster_rows = F,cluster_cols = F,show_rownames = T,show_colnames = F,border_color = "white",color = colorRampPalette(colors = c("#5EAD83","white","#9F558E"))(100))
##########################################################













