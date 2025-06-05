rm(list=ls())
setwd('/Users/dichen/Documents/Methylation_Running/')

library(R.matlab) 
library('CCA')
library(foreach)
library(doParallel)
library(ggExtra)

cd_Eta_adj <- function(Eta_real,Eta_10000){
  Eta_0<-mean(Eta_10000)
  Eta_adj<-1-(1-Eta_real)/(1-Eta_0)
  return (Eta_adj)
}

ES_CC <- function(Cor) {
  Phi <- Cor^2
  Lambda <- exp(sum(log(1-Phi)))
  Eta <- 1-Lambda
  HT <- sum(Phi/(1-Phi))
  Output <- list(Eta=Eta, HT=HT)
  return (Output)
}

PermFun <- function(M1,M2,Regular){
  library('CCA')
  Perm <- sample(nrow(M2),replace = F)
  Perm_RCC <- rcc(M1,M2[Perm,],Regular,Regular)
  Perm_Cor <- Perm_RCC$cor
  Perm_Eta <- (ES_CC(Perm_Cor))$Eta
  return (c(Perm_Eta,Perm_Cor))
}

######## load data
select_cluster<-c(2,10,13:18);# select brain 8 clusters

length(select_cluster)
fMRI_type<- c("SurfArea","ThickAvg","SubcortexVol")

for (k in 1:3) {
# k<-1
print(fMRI_type[k])

data_dir<-paste0("./data_CCA_change_Methy_",fMRI_type[k],".mat")
data<-readMat(data_dir, header = TRUE)
data<-readMat(data_dir)

ME_change<-data$ME.change.regressed
ME_change_select_cluster9<-ME_change[,select_cluster]# selcet 10 clusters

T1_change<-data$T1.change.regressed
sample_size<-nrow(T1_change)

############# parameter
Regular <-0 ## ridges
NPerm <- 10000 ## permutation times

########### CCA
RCC <- rcc(ME_change_select_cluster9,T1_change,Regular,Regular)
Origin_Cor <- RCC$cor
Origin_Eta <- (ES_CC(RCC$cor))$Eta
#assign ('Eta', Eta)
cor(RCC$scores$xscores[,1],RCC$scores$yscores[,1])

RCC$scores$corr.X.xscores
RCC$scores$corr.Y.yscores

RCC$scores$xscores[1:3,1:3]
kkk<-ME_change_select_cluster9 %*% RCC$xcoef

cor(RCC$scores$xscores[1,],kkk[1,])

library(psych)

## r-matrix ME
r_ME_all<-matrix(data = NA, nrow = 3, ncol = 8)
p_ME_all<-matrix(data = NA, nrow = 3, ncol = 8)
for (i in 1:3){
  i
for (j in 1:8){
resul_cor<-corr.test(RCC$scores$xscores[,i],ME_change_select_cluster9[,j],method = "pearson") 
p_ME_all[i,j]<-resul_cor$p
r_ME_all[i,j]<-resul_cor$r
}
}
ME_r<-t(r_ME_all)
ME_p<-t(p_ME_all)

## r-matrix MRI
r_T1_all<-matrix(data = NA, nrow = 3, ncol = ncol(T1_change))
p_T1_all<-matrix(data = NA, nrow = 3, ncol = ncol(T1_change))
for (i in 1:3){
  i
  for (j in 1:ncol(T1_change)){
    #r<-cor(RCC$scores$yscores[,i],T1_change[,j])
    
    resul_cor<-corr.test(RCC$scores$yscores[,i],T1_change[,j],method = "pearson") 
    p_T1_all[i,j]<-resul_cor$p
    r_T1_all[i,j]<-resul_cor$r
  }
}
MRI_r<-t(r_T1_all)
MRI_p<-t(p_T1_all)

############### Permutation
Perm_output <- foreach(x=1:NPerm,.combine='rbind') %dopar% PermFun(ME_change_select_cluster9,T1_change,Regular)
#dim(Perm_output)
P_perm <- rowMeans(apply(Perm_output,1,">",c(Origin_Eta,Origin_Cor)))

P_perm
Origin_Eta
cd_Eta_adj(Origin_Eta,Perm_output[,1])
############### load saved p_permu

print(fMRI_type[k])
result_dir<-paste0("./Result_CCA_Nobrain_cluster8/result_CCA_change_Methy_",fMRI_type[k],".mat")
output_load<-readMat(result_dir)

output_load$P.perm
output_load$Origin.Eta
p.adjust(output_load$P.perm[2:4], method = "BH")
sum(output_load$Perm.output[,1]>Origin_Eta)/10000

## adjust Eta
cd_Eta_adj(output_load$Origin.Eta,output_load$Perm.output[,1])

##################### PLot corr
i<-1# ith component

x<-RCC$scores$xscores[,i]
y<-RCC$scores$yscores[,i]
cor(x,y)

## dele outlier
#xy<-cbind(x,y)
#xy<-xy[xy[,2] <= 5,]
#xy<-xy[xy[,2] >-3,]

#x<-xy[,1]
#y<-xy[,2]
#cor(x,y)

##plot
data1 <- data.frame(x=x,y=y)
colors_1 <- c("#e3dae6", "#B3CB97","#7995A2")
#(a) 二维散点与统计直方图
scatter <- ggplot(data=data1,aes(x=x,y=y)) + 
  geom_point(shape=21,fill=colors_1[k],color=colors_1[k],size=3,alpha=0.66)+
  theme_minimal()+
  geom_smooth(formula = 'y ~ x',method="lm", color="black", size=1.5) +  # 添加拟合直线
  theme(
    #text=element_text(size=15,face="plain",color="black"),
    #axis.title=element_text(size=15,face="plain",color="black"),
    axis.title.x=element_blank(),  # 移除 x 轴标题
    axis.title.y=element_blank(),  # 移除 x 轴标题
    axis.text = element_text(size=13,face="plain",color="black"),
    legend.text= element_text(size=13,face="plain",color="black"),
    legend.title=element_text(size=12,face="plain",color="black"),
    axis.line = element_line(color="black", size=0.5),  # 设置轴线颜色和大小
    legend.background=element_blank()
    #legend.position = c(0.12,0.88)
  )
ggMarginal(scatter,type="density",color="black",fill=colors_1[k])
print(fMRI_type[k])
##################### PLot density
library(ggplot2)
Eta_permution<-output_load$Perm.output[,1]
Eta_permution<-data.frame(Eta_permution)
colnames(Eta_permution)<-c("Eta_10000")

## define color
colors_1 <- c("#e3dae6", "#B3CB97","#7995A2")
colors_2 <- c("#8b80a8", "#6B7E53","#284C69")

## plot
ggplot(Eta_permution, aes(x=Eta_10000))+ 
  geom_histogram(aes(y=..density..), colour=colors_2[k],fill=colors_1[k])+
  geom_density(fill = colors_1[k], alpha = 0.7, color = colors_2[k])+
  #geom_histogram(alpha=0.55,bw=1,colour="black",size=0.25)+
  
  labs(x = "Eta-squared value", y = "Density") +
  #theme_classic()+
  theme(plot.title = element_text(size = 20, face = "bold"),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 16),
        panel.background = element_rect(fill = "#f7f7f7"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  geom_vline(xintercept = output_load$Origin.Eta, linetype = "solid", color = "#d54846",size = 0.8)
print(fMRI_type[k])

output_load$Origin.Eta
output_load$P.perm[1]
##################### save component
save_component_dir<-paste0("./Result_CCA_Nobrain_cluster8/Component_",fMRI_type[k],".mat")
ME_Component_123<-RCC$scores$xscores[,1:3]
T1_Component_123<-RCC$scores$yscores[,1:3]
#writeMat(save_component_dir,ME_Component_123=ME_Component_123,T1_Component_123=T1_Component_123)

##################### save
save_dir<-paste0("./Result_CCA_Nobrain_cluster8/result_CCA_change_Methy_",fMRI_type[k],".mat")
#writeMat(save_dir,Perm_output=Perm_output, sample_size=sample_size,P_perm=P_perm,Origin_Eta=Origin_Eta,Origin_Cor=Origin_Cor)

}

##################### 
## adj the three Eta P-value
k=1
result_dir<-paste0("./Result_CCA_Nobrain_cluster8/result_CCA_change_Methy_",fMRI_type[k],".mat")
output_load<-readMat(result_dir)
P_1<-output_load$P.perm[1]

k=2
result_dir<-paste0("./Result_CCA_Nobrain_cluster8/result_CCA_change_Methy_",fMRI_type[k],".mat")
output_load<-readMat(result_dir)
P_2<-output_load$P.perm[1]

k=3
result_dir<-paste0("./Result_CCA_Nobrain_cluster8/result_CCA_change_Methy_",fMRI_type[k],".mat")
output_load<-readMat(result_dir)
P_3<-output_load$P.perm[1]

c(P_1,P_2,P_3)

## FDR
P_123<-c(P_1,P_2,P_3)
p.adjust(P_123, method = "BH")
##################### 
# finished
################ 












