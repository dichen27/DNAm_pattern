
rm(list=ls())
setwd('/Users/dichen/Documents/Methylation_Running/')

library(vegan)
library(R.matlab)
library(pheatmap)
##########################
## load IMAGEN regress pattern (final pattern)
mat_data <- readMat("./(Reply1-8)pattern_longitudinal_cluster18_regressed.mat")
r_IMAGEN_change<- mat_data$pattern.longitudinal.cluster18.regressed

## heatmap-change-cluster18
data_change<- readMat("./(Reply1-8)substract_cluster18_regress_19_14.mat")
data_cluster18_change_regress<-data_change$substract.cluster18.regress.19.14

## order the cluster
num_order<-c(4,1,13,12,6,11,18,8,5,2,14,3,15,10,9,7,16,17)

data_change_cluster18_order<-matrix(data=NA, nrow =506, ncol =18, byrow = FALSE, dimnames = NULL)
for (i in 1:18){
  k<-num_order[i]
  data_change_cluster18_order[,i]<-data_cluster18_change_regress[,k]
}

r_change_cluster18_order<- cor(data_change_cluster18_order, method = "pearson", use = "complete.obs") 

rnames<-paste0("Cluster-", num_order)
cnames<-paste0("Cluster-", num_order)
r_data<-matrix(r_change_cluster18_order, nrow = 18,ncol = 18,dimnames=list(rnames, cnames))
pheatmap(r_data, cluster_rows = F,cluster_cols = F,show_rownames = T,show_colnames = F,border_color = "white",color = colorRampPalette(colors = c("#5EAD83","white","#9F558E"))(100))

##########################
## load pattern 14
pattern_14 <- readMat("./(Reply1-8)pattern_age14_regress.mat")
r_pattern_14_regress<- pattern_14$pattern.age14.regress

# Mantel test
mantel(xdis = r_IMAGEN_change, ydis = r_pattern_14_regress, method = "pearson", permutations = 10000, na.rm = TRUE)

## heatmap-14
data_14<- readMat("./(Reply1-8)data_cluster18_age14_regressed.mat")
data_cluster18_age14_regress<-data_14$data.cluster18.age14.regressed

## order the cluster
num_order<-c(4,1,13,12,6,11,18,8,5,2,14,3,15,10,9,7,16,17)

data_14_cluster18_order<-matrix(data=NA, nrow =506, ncol =18, byrow = FALSE, dimnames = NULL)
for (i in 1:18){
  k<-num_order[i]
  data_14_cluster18_order[,i]<-data_cluster18_age14_regress[,k]
}

r_14_cluster18_order<- cor(data_14_cluster18_order, method = "pearson", use = "complete.obs") 

rnames<-paste0("Cluster-", num_order)
cnames<-paste0("Cluster-", num_order)
r_data<-matrix(r_14_cluster18_order, nrow = 18,ncol = 18,dimnames=list(rnames, cnames))
pheatmap(r_data, cluster_rows = F,cluster_cols = F,show_rownames = T,show_colnames = F,border_color = "white",color = colorRampPalette(colors = c("#5EAD83","white","#9F558E"))(100))


##########################
## load pattern 19
pattern_19 <- readMat("./(Reply1-8)pattern_age19_regress.mat")
r_pattern_19_regress<- pattern_19$pattern.age19.regress

# Mantel test
mantel(xdis = r_IMAGEN_change, ydis = r_pattern_19_regress, method = "pearson", permutations = 10000, na.rm = TRUE)

############################ heatmap-19
data_19<- readMat("./(Reply1-8)data_cluster18_age19_regressed.mat")
data_cluster18_age19_regress<-data_19$data.cluster18.age19.regressed

## order the cluster
num_order<-c(4,1,13,12,6,11,18,8,5,2,14,3,15,10,9,7,16,17)

data_19_cluster18_order<-matrix(data=NA, nrow =506, ncol =18, byrow = FALSE, dimnames = NULL)
for (i in 1:18){
  k<-num_order[i]
  data_19_cluster18_order[,i]<-data_cluster18_age19_regress[,k]
}

r_19_cluster18_order<- cor(data_19_cluster18_order, method = "pearson", use = "complete.obs") 

r_data<-matrix(r_19_cluster18_order, nrow = 18,ncol = 18,dimnames=list(rnames, cnames))
pheatmap(r_data, cluster_rows = F,cluster_cols = F,show_rownames = T,show_colnames = F,border_color = "white",color = colorRampPalette(colors = c("#5EAD83","white","#9F558E"))(100))

##########################
## finished
