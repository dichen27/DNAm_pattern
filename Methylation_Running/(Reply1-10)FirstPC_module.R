
rm(list=ls())
setwd('/Users/dichen/Documents/Methylation_Running/')
#setwd('F:/Methylation_Running/')
library(corpcor)
library(R.matlab)
########################################################################
num <- 372582
sub <- 506
mat_data <- readMat("./List_ME/substract_19_14_all_CpGs.mat")
substract_19_14_all_CpGs<-mat_data$substract.19.14.all.CpGs

data_substract<-substract_19_14_all_CpGs[,-1]

########################################################################
load('./IMAGEN_methylation_data_origin/t14_506_ori.rda')
all_cpg_IMAGEN<-colnames(t14_2)
all_cpg_IMAGEN<-all_cpg_IMAGEN[2:length(all_cpg_IMAGEN)] # dele the "ID"


## load ME175
bwnet <- readRDS('./Consensus_19_14_Network_Information_Power5.rds')
Origin <- as.data.frame(bwnet$colors)
index_CpG<-Origin$`bwnet$colors`


length(unique(index_CpG))
sum(index_CpG==3)
sig_cpg_cmd<-all_cpg_IMAGEN[index_CpG==3]
length(sig_cpg_cmd)

########################################################################
## PCA
Variance_explained_first_PC<-matrix(data=NA, nrow = 175, ncol =1, byrow = FALSE, dimnames = NULL)
num_CpGs_in_module<-matrix(data=NA, nrow = 175, ncol =1, byrow = FALSE, dimnames = NULL)
for (k in 1:175){
  print(k)
  
  ID_CpGs_k<-all_cpg_IMAGEN[index_CpG==k]
  data_k<-matrix(data=NA, nrow = 506, ncol =length(ID_CpGs_k), byrow = FALSE, dimnames = NULL)
  for (i in 1:length(ID_CpGs_k)){
    lin<-data_substract[,all_cpg_IMAGEN==ID_CpGs_k[i]]
    data_k[,i]<-lin
  }
  
  ##PCA
  ss <- fast.svd(data_k,tol=0) # this may take a while
  percvar <- ss$d^2/sum(ss$d^2)
  # plot(percvar,xlab="Principal Components",ylab="Variance explained")
 
  Variance_explained_first_PC[k]<-percvar[1]
  num_CpGs_in_module[k]<-length(ID_CpGs_k)
}
writeMat("./(Reply1-10)Variance_explained_first_PC.mat", Variance_explained_first_PC = Variance_explained_first_PC)
writeMat("./(Reply1-10)num_CpGs_in_module.mat", num_CpGs_in_module = num_CpGs_in_module)

########################################################################
## analysis result
result<-readMat("./(Reply1-10)Variance_explained_first_PC.mat")
Variance_explained_first_PC<-result$Variance.explained.first.PC

n<-175
mean_data<-mean(Variance_explained_first_PC)
mean_data

sd_data<-sd(Variance_explained_first_PC)
sd_data

mean_data-(1.96*sd_data)/sqrt(n)
mean_data+(1.96*sd_data)/sqrt(n)


## 95% CI
t_critical <- qt(0.975, df = n - 1)  # two-tailed, hence 0.975

# Calculate the margin of error
margin_error <- t_critical * sd_data / sqrt(n)

# Calculate the confidence interval
ci_lower <- mean_data - margin_error
ci_lower

ci_upper <- mean_data + margin_error
ci_upper



########################################################################
## finished
