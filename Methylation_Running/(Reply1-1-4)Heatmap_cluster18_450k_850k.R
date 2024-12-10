
rm(list=ls())
 setwd('/Users/dichen/Documents/Methylation_Running/')


library(R.matlab)
library(flashClust)
library(pheatmap)
library(vegan)
 
 ## load ME_module190_substract_19_14 n=506
ME_module175_substract_19_14 <- read.csv('./List_ME/ME_consensus_substract_19_14.csv')
 

r_substract <- read.csv('./List_ME/r_substract.csv')      
sampleTree <- flashClust(dist(r_substract,method = 'euclidean'), method="average")
order_sampleTree<-as.double(sampleTree$order) #num of ME, used for every heatmap

###########################################################################
# Heatmap change n=506
ME_module175<-ME_module175_substract_19_14
IMAGEN_weight_order<-matrix(data=NA, nrow =506, ncol =175, byrow = FALSE, dimnames = NULL)
for (i in 1:175){

  k<-order_sampleTree[i]
  IMAGEN_weight_order[,i]<-ME_module175[,k]
}
r_substract_19_14_order<- cor(IMAGEN_weight_order, method = "pearson", use = "complete.obs") 
pheatmap(r_substract_19_14_order, cluster_rows = F,cluster_cols = F,show_rownames = F,show_colnames = F,border_color = "white",color = colorRampPalette(colors = c("#5EAD83","white","#9F558E"))(100))

###########################################################################
# Heatmap change 450k
ME_module175<-read.csv('./List_ME/ME_consensus_substract_19_14_450k.csv')
IMAGEN_weight_order<-matrix(data=NA, nrow =445, ncol =175, byrow = FALSE, dimnames = NULL)
for (i in 1:175){
  
  k<-order_sampleTree[i]
  IMAGEN_weight_order[,i]<-ME_module175[,k]
}
r_change_450k<- cor(IMAGEN_weight_order, method = "pearson", use = "complete.obs") 
pheatmap(r_change_450k, cluster_rows = F,cluster_cols = F,show_rownames = F,show_colnames = F,border_color = "white",color = colorRampPalette(colors = c("#5EAD83","white","#9F558E"))(100))
# mantel
kkk<-mantel(xdis = r_change_450k, ydis = r_substract_19_14_order, method = "pearson", permutations = 10000, na.rm = TRUE)
kkk$statistic


###########################################################################
# Heatmap change 850k
ME_module175<-read.csv('./List_ME/ME_consensus_substract_19_14_850k.csv')
IMAGEN_weight_order<-matrix(data=NA, nrow =61, ncol =175, byrow = FALSE, dimnames = NULL)
for (i in 1:175){
  
  k<-order_sampleTree[i]
  IMAGEN_weight_order[,i]<-ME_module175[,k]
}
r_change_850k<- cor(IMAGEN_weight_order, method = "pearson", use = "complete.obs") 
pheatmap(r_change_850k, cluster_rows = F,cluster_cols = F,show_rownames = F,show_colnames = F,border_color = "white",color = colorRampPalette(colors = c("#5EAD83","white","#9F558E"))(100))

# mantel
kkk<-mantel(xdis = r_change_850k, ydis = r_substract_19_14_order, method = "pearson", permutations = 10000, na.rm = TRUE)
kkk$statistic

# mantel with 450k
kkk<-mantel(xdis = r_change_850k, ydis = r_change_450k, method = "pearson", permutations = 10000, na.rm = TRUE)
kkk$statistic










