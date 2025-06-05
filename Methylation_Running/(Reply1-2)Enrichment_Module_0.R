rm(list=ls())
setwd('/Users/dichen/Documents/Methylation_Running/')

library(missMethyl)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(limma)
library(ggplot2)
library(pheatmap)
library("TissueEnrich")
library("FDb.InfiniumMethylation.hg19")
hm450.hg19 <- get450k(genome='hg19')
###################### missMehtyl example data
ann <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)

# All CpG sites tested
allcpgs <- rownames(ann)
allcpgs[1]

##################### load all IMAGEN CpGs
# age 14
load('./IMAGEN_methylation_data_origin/t14_506_ori.rda')
t14_2[1:3,1:3]
ID_506<-matrix(t14_2[,1])


all_cpg_IMAGEN<-colnames(t14_2)
#all_cpg_IMAGEN<-colnames(all_cpg_IMAGEN)

all_cpg_IMAGEN[1:3]
all_cpg_IMAGEN<-all_cpg_IMAGEN[2:length(all_cpg_IMAGEN)] # dele the "ID"

length(all_cpg_IMAGEN)
all_cpg_IMAGEN[1:3]

################### load the Module-0

list_all <- read.table("./List_Of_Methylation_Change_Power5_Min50_Merge0_UnsignedTOM_Consensus_19_14.txt", sep = "\t")  # Adjust sep if necessary
length(list_all$V1)
list_all<-list_all$V1
list_all[1:6]
tabulate(list_all)

index_module_0<-list_all==0;

## Number of Module_0 CpGs
sum(index_module_0)

CpGs_module_0<-all_cpg_IMAGEN[index_module_0]

saveRDS(CpGs_module_0, file="./(Reply1-2)CpGs_module_0.rds",compress=T)

########################
## data_module_0 at age-14

t14_matrix<-t14_2[,2:ncol(t14_2)]
dim(t14_matrix)
t14_matrix[1:3,1:3]

data_14_module_0<-matrix(data=NA, nrow = 506, ncol =1, byrow = FALSE, dimnames = NULL)
ID_CpGs_k <- readRDS("./(Reply1-2)CpGs_module_0.rds")
  
  data_k<-matrix(data=NA, nrow = 506, ncol =length(ID_CpGs_k), byrow = FALSE, dimnames = NULL)
  for (i in 1:length(ID_CpGs_k)){
    lin<-t14_matrix[,colnames(t14_matrix)==ID_CpGs_k[i]]
    data_k[,i]<-lin
  }
  
data_14_CpG_based_module_0<-data_k

save(data_14_CpG_based_module_0, file="./(Reply1-2)data_14_CpG_based_module_0.rda")
  
data_k_mean<-apply(data_k,1,mean)  # mean of the each row
data_14_module_0<-data_k_mean
write.csv(data_14_module_0, file = "./(Reply1-2)data_14_module_0.csv")
###################### age-19
load('./IMAGEN_methylation_data_origin/t19_506_ori.rda')
t19_2[1:3,1:3]
dim(t19_2)
t19_matrix<-t19_2[,2:ncol(t19_2)]
dim(t19_matrix)
t19_matrix[1:3,1:3]


data_19_module_0<-matrix(data=NA, nrow = 506, ncol =1, byrow = FALSE, dimnames = NULL)
ID_CpGs_k <- readRDS("./(Reply1-2)CpGs_module_0.rds")
  
  data_k<-matrix(data=NA, nrow = 506, ncol =length(ID_CpGs_k), byrow = FALSE, dimnames = NULL)
  for (i in 1:length(ID_CpGs_k)){
    lin<-t19_matrix[,colnames(t19_matrix)==ID_CpGs_k[i]]
    data_k[,i]<-lin
  }

data_19_CpG_based_module_0<-data_k
save(data_19_CpG_based_module_0, file="./(Reply1-2)data_19_CpG_based_module_0.rda")


  data_k_mean<-apply(data_k,1,mean)  # mean of the each row
  data_19_module_0<-data_k_mean

write.csv(data_19_module_0, file = "./(Reply1-2)data_19_module_0.csv")

#################### Enrichment
########################
# GO testing with prior probabilities taken into account
gst <- gometh(sig.cpg = CpGs_module_0, all.cpg = all_cpg_IMAGEN, collection = "GO", 
              plot.bias = TRUE, prior.prob = TRUE, anno = ann)
# Total number of GO categories significant at 5% FDR
table(gst$FDR<0.05)
# Table of top GO results
topGSA(gst,n=10)


#### save GO
gst_dir<-"./Enrichment_Second_Round_Cluster18/Enrichment_GO_Module_0.rds"
saveRDS(gst, file=gst_dir,compress=T)

########## plot GO
## read
go <- readRDS(gst_dir)
topGSA(go,n=10)

top_10<-topGSA(go,n=10)
index<-top_10$FDR<0.05
top_go<-top_10[index,]

## plot
Description <-top_go$TERM
P_FDR <-top_go$FDR
Gene_Count <-top_go$N
df <- data.frame(Description, P_FDR, Gene_Count)

ggplot(df, aes(x = -log10(P_FDR), y = Description, size = Gene_Count)) +
  #scale_fill_manual(values=c(c = "red", d = "blue", e = "green" , p = "orange", r = "yellow"))+
  geom_point(alpha=0.97,color = "#7670AA") +
  scale_size(range = c(2, 6)) +
  labs(x = "-Log10(FDR P-value)", y = " ", title = "") +
  theme_minimal()+
  theme(plot.margin = margin(10, 10, 10, 10, "pt"),
        panel.border = element_rect(color = "gray77", fill = NA, linewidth = 1.2),
        axis.text.x = element_text(size = 10, color = "black"),
        axis.text.y = element_text(size = 10, color = "black"))

########################################################################
# KEGG testing
kegg <- gometh(sig.cpg = CpGs_module_0, all.cpg = all_cpg_IMAGEN, collection = "KEGG", 
               prior.prob=TRUE, anno = ann)
# FDR
table(kegg$FDR<0.05)
# Table of top KEGG results
topGSA(kegg,n=10)


#### save KEGG
kegg_dir<-"./Enrichment_Second_Round_Cluster18/Enrichment_KEGG_Module_0.rds"
saveRDS(kegg, file=kegg_dir,compress=T)
########################################################################
########################################################################
## finished










