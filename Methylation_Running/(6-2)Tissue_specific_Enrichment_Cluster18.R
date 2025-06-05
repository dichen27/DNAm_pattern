### for tissue specific gene enrichment analysis
### for the first time to run this analysis, install them via the following arguments
rm(list=ls())
setwd('/Users/dichen/Documents/Methylation_Running/')


### to load the package
library("TissueEnrich")
library("FDb.InfiniumMethylation.hg19")
hm450.hg19 <- get450k(genome='hg19')

for (k in 1:18) {
#k<-9
print(k)
### to load the probe names
cluster5_dir<-paste(getwd(),"/Cluster18/CLuster_CpGs_k",k,".rds",sep = "")

sig_cpg_IMAGEN_k<- readRDS(cluster5_dir)
length(sig_cpg_IMAGEN_k)

### to load the probe names
moduleProbes <- sig_cpg_IMAGEN_k

### to get their genomic locations
anno_moduleProbes <- hm450.hg19[c(moduleProbes)]

### to get their nearest genes
TSS_moduleProbes <- getNearestTSS(anno_moduleProbes)

### to remove the overlapping records
Gene_module <- unique(TSS_moduleProbes$nearestGeneSymbol)

### to perform the tissue-specific enrichment analysis
gs<-GeneSet(geneIds=Gene_module,organism='Homo Sapiens',
            geneIdType=SymbolIdentifier())
output<-teEnrichment(gs)
seEnrichmentOutput<-output[[1]]
enrichmentOutput<-setNames(data.frame(assay(seEnrichmentOutput),
                                      row.names = rowData(seEnrichmentOutput)[,1]),
                           colData(seEnrichmentOutput)[,1])
enrichmentOutput$Tissue<-row.names(enrichmentOutput)

### to save the results
save_dir<-paste(getwd(),"/Enrichment_Tissue_specific_result_Cluster18/Enrichment_Tissue_specific_k",k,".RData",sep = "")
save(enrichmentOutput, file = save_dir)
}

############################################################################################################################
### to plot the results
# load Enrichment-result
k<-12
dir<-paste(getwd(),"/Enrichment_Tissue_specific_result_Cluster18/Enrichment_Tissue_specific_k",k,".RData",sep = "")
load(dir)

# plot Enrichment-result
title_name<-paste("Tissue-specific enrichment for Cluster-",k,sep = "")

ggplot(enrichmentOutput,aes(x=reorder(Tissue,-Log10PValue),y=Log10PValue,
                            label = Tissue.Specific.Genes,fill = Tissue))+
  geom_bar(stat = 'identity')+
  labs(x="", y = '-LOG10(P-Value)')+
  ggtitle(title_name)+
  theme_bw()+
  theme(legend.position='none')+
  theme(plot.title = element_text(hjust = 0.5,size = 15),axis.title =
          element_text(size=10))+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        panel.grid.major= element_blank(),panel.grid.minor = element_blank())+
  theme(plot.margin = unit(c(1,1,1,1), "cm"))

## statistics of Cerebral Cortex
index<-enrichmentOutput$Tissue=="Cerebral Cortex"
P_log<-enrichmentOutput$Log10PValue[index]

## Bonferroni corrected P-value
(10^-P_log)*length(enrichmentOutput$Tissue)

#####################################################################
## write the Enrichment-function
cd_Enrichment_Tissue_specific <- function(sig_cpg_IMAGEN){
  library("TissueEnrich")
  library("FDb.InfiniumMethylation.hg19")
  hm450.hg19 <- get450k(genome='hg19')
  
  ### to load the probe names
  moduleProbes <- sig_cpg_IMAGEN
  
  ### to get their genomic locations
  anno_moduleProbes <- hm450.hg19[c(moduleProbes)]
  
  ### to get their nearest genes
  TSS_moduleProbes <- getNearestTSS(anno_moduleProbes)
  
  ### to remove the overlapping records
  Gene_module <- unique(TSS_moduleProbes$nearestGeneSymbol)
  
  ### to perform the tissue-specific enrichment analysis
  gs<-GeneSet(geneIds=Gene_module,organism='Homo Sapiens',
              geneIdType=SymbolIdentifier())
  output<-teEnrichment(gs)
  seEnrichmentOutput<-output[[1]]
  enrichmentOutput<-setNames(data.frame(assay(seEnrichmentOutput),
                                        row.names = rowData(seEnrichmentOutput)[,1]),
                             colData(seEnrichmentOutput)[,1])
  enrichmentOutput$Tissue<-row.names(enrichmentOutput)
  
  return (enrichmentOutput)
}

## try the new function
enrichmentOutput<-cd_Enrichment_Tissue_specific(sig_cpg_IMAGEN_k)



