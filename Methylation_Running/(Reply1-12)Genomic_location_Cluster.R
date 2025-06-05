### for tissue specific gene enrichment analysis
### for the first time to run this analysis, install them via the following arguments
rm(list=ls())
setwd('/Users/dichen/Documents/Methylation_Running/')


### to load the package
library("FDb.InfiniumMethylation.hg19")
hm450.hg19 <- get450k(genome='hg19')


Cluster_plot_all<-matrix(NA, nrow = 1, ncol = 4)
CHR_percent_all_cluster<-matrix(NA, nrow = 1, ncol = 22)
## find the nearest Gene
for (k in 1:18) {

# k<-1
print(k)
### to load the probe names
cluster5_dir<-paste(getwd(),"/Cluster18/CLuster_CpGs_k",k,".rds",sep = "")

sig_cpg_IMAGEN_k<- readRDS(cluster5_dir)
print(length(sig_cpg_IMAGEN_k))

### to load the probe names
moduleProbes <- sig_cpg_IMAGEN_k

### to get their genomic locations
anno_moduleProbes <- hm450.hg19[c(moduleProbes)]
#anno_moduleProbes <- hm450.hg19["ch.9.22352682F"]

CHR<-as.character(anno_moduleProbes)

## map to plot document
CHR_num <- sub("chr(\\d+):.*", "\\1", CHR)

POS <- sub(".*:(\\d+)-.*", "\\1", CHR)
POS_cleaned <- sub(".*:", "", POS)


PHENOTYPE<-rep(paste0("Cluster-",k), length(sig_cpg_IMAGEN_k))
CpG_ID<-names(CHR)
GROUPCOLOR<-rep(paste0("Cluster-",k), length(sig_cpg_IMAGEN_k))

#index<-CHR_num==1
#sum(index)

Cluster_plot_k<-cbind(CHR_num,POS_cleaned,PHENOTYPE,GROUPCOLOR)
colnames(Cluster_plot_k) <- c("CHR","POS","PHENOTYPE","GROUPCOLOR")


set.seed(123) 
if (nrow( Cluster_plot_k) >= 10000) {
  Cluster_plot_k_select <-Cluster_plot_k[sample(nrow(Cluster_plot_k), 10000), ]
} else {
  Cluster_plot_k_select<- Cluster_plot_k
}

#set.seed(123)
#sampled_Cluster_plot_k <- Cluster_plot_k[sample(nrow(Cluster_plot_k), size = 0.2 * nrow(Cluster_plot_k)), ]


### save Cluster-k
save_dir<-paste0("./Cluster_plot/Cluster_plot_k",k,".txt")
write.table(Cluster_plot_k_select, file = save_dir, sep = "\t", row.names = FALSE)

## CHR_percent
freq_table <- table(CHR_num)
percent_table <- prop.table(freq_table)
percent_sorted <- percent_table[order(as.numeric(names(percent_table)))]
CHR_percent_all_cluster<-rbind(CHR_percent_all_cluster,percent_sorted)


### all Cluster
Cluster_plot_all<-rbind(Cluster_plot_all,Cluster_plot_k)
}


##
CHR_percent_all_cluster<-CHR_percent_all_cluster[-1, ] # dele row-1

##########################################################################################
## finished
