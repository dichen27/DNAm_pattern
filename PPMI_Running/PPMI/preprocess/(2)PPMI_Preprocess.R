
rm(list=ls())
library(minfi)
library(R.matlab) 
setwd("/Users/dichen/Documents/PPMI/preprocess/")

load("./RGset_PPMI.rda")
pd <- pData(RGset)
qcReport(RGset, sampNames = pd$rownames, sampGroups = NULL, pdf = "./figures/qcReport_PPMI_No_normalize.pdf")

pdf(file="./figures/BS_conversion_No_normalize.pdf")
	controlStripPlot(RGset, controls="BISULFITE CONVERSION II", sampNames = pd$IMAGEN_PSC1)
dev.off()



object=preprocessIllumina(RGset) # preprocess
object<-mapToGenome(object)
object=ratioConvert(object, type="Illumina") #convert raw methylation data
beta <- getBeta(object) #get the beta value for each probe
dat <- object
pd=pData(dat)

dim(beta)
beta[1:3,1:3]


##################################
## Use MDS to check bad samples ##
##################################

###### load self-report gender_order
data_self<-readMat("gender_self_order.mat", header = TRUE)
gender_self_order<-data_self$gender.self.order
gender_self<-gender_self_order[,2]
tabulate(gender_self)


keepIndex=which(!seqnames(dat)%in%c("chrX","chrY"))
beta <- getBeta(dat)
beta <- beta[keepIndex,]
mydist=dist(t(beta[sample(1:nrow(beta),10000),]))
mds=cmdscale(mydist)
length(mds[,1])
pdf("./figures/gender_effect_mds.pdf",width=25,height=4)
	# plot(mds[,1],xlab="",ylab="First Principal component",xaxt="n")
plot(mds[,1],col=as.numeric(as.factor(gender_self))+1,xlab="",ylab="First Principal component",xaxt="n")
dev.off()

###### Race effect
data_self<-readMat("Race_self_order.mat", header = TRUE)
Race_self_order<-data_self$Race.self.order
Race_self<-Race_self_order[,2]
length(Race_self)


keepIndex=which(!seqnames(dat)%in%c("chrX","chrY"))
beta <- getBeta(dat)
beta <- beta[keepIndex,]
mydist=dist(t(beta[sample(1:nrow(beta),10000),]))
mds=cmdscale(mydist)
length(mds[,1])
pdf("./figures/Race_effect_mds.pdf",width=25,height=4)
# plot(mds[,1],xlab="",ylab="First Principal component",xaxt="n")
plot(mds[,1],col=as.numeric(as.factor(Race_self))+1,xlab="",ylab="First Principal component",xaxt="n")
dev.off()

####################################
## Use PCA to check batch effects ##
####################################
b <- beta - rowMeans(beta)
# Do the svd 
library(corpcor)
fileName <-"./fast_svd.rda"
if (file.exists(fileName)) {
    load(fileName)
} else {
	ss <- fast.svd(b,tol=0) # this may take a while
	save(ss,file="./fast_svd.rda")
}
# Look at the percent of variance explained
# by each principal component
percvar <- ss$d^2/sum(ss$d^2)
pdf(file="./figures/PCA_distribution.pdf")
	plot(percvar,xlab="Principal Components",ylab="Variance explained")
dev.off()

################################### PCs_gender_effect
	pdf(file=paste("./figures/PCs_gender_effect.pdf",sep=""))
	variabletocolor=gender_self
	bob=levels(factor(variabletocolor))
	colors=match(variabletocolor,bob)
	pairs(ss$v[,1:4], col=colors, labels=c("PC1", "PC2", "PC3", "PC4"))
	par(mfrow=c(2,2))
	boxplot(ss$v[,1]~gender_self,ylab="PC1",xlab="Gender")
	boxplot(ss$v[,2]~gender_self,ylab="PC2",xlab="Gender")
	boxplot(ss$v[,3]~gender_self,ylab="PC3",xlab="Gender")
	boxplot(ss$v[,4]~gender_self,ylab="PC4",xlab="Gender")
	dev.off()


	################################### PCs_slide_effect
	ID_524<-object@colData@rownames
	ID_slide<-substr(ID_524,1,12)
	
	
	pdf(file=paste("./figures/PCs_slide_effect.pdf",sep=""))
	variabletocolor=ID_slide
	bob=levels(factor(variabletocolor))
	colors=match(variabletocolor,bob)
	pairs(ss$v[,1:4], col=colors, labels=c("PC1", "PC2", "PC3", "PC4"))
	par(mfrow=c(2,2))
	boxplot(ss$v[,1]~ID_slide,ylab="PC1",xlab="Slide")
	boxplot(ss$v[,2]~ID_slide,ylab="PC2",xlab="Slide")
	boxplot(ss$v[,3]~ID_slide,ylab="PC3",xlab="Slide")
	boxplot(ss$v[,4]~ID_slide,ylab="PC4",xlab="Slide")
	dev.off()
	
	### Mark individuals out of normal range (i.e. median+3SD or median-3SD) based on first 4 components ###
	RM <- rep(FALSE, nrow(ss$v))
	for (i in 1:4) {
	  Median <- median(ss$v[,i]) #calculate median of each component
	  SD <- sd(ss$v[,i]) #calculate the standard deviation of each component
	  RM <- RM|(abs(ss$v[,i]-Median)>4*SD) #mark individuals outside 3SD range as TRUE
	}
	# RM will be used to remove outliers at a later stage
sum(RM)



############ predicted_gender
predicted_gender <- getSex(object)
predicted_gender <- as.numeric(as.factor(predicted_gender[,3]))

object@rowRanges[233]
ID_524<-object@colData@rownames
writeMat("./ID_524.mat",ID_524=ID_524)
tabulate(predicted_gender)

## load self-report gender_order
data_self<-readMat("gender_self_order.mat", header = TRUE)
gender_self_order<-data_self$gender.self.order
gender_self<-gender_self_order[,2]
tabulate(gender_self)


############ Identify samples with gender discrepancy
pdf(file="./figures/Gender_QC.pdf")
plot(jitter(as.numeric(factor(gender_self)))~jitter(predicted_gender), xlab="Gender (predicted)", ylab="Gender (Reported)", main="QC (gender discrepancy)",xaxt="n",yaxt="n")
axis(2,c(1,2),c("Female","Male"))
axis(1,c(1,2),c("Female","Male"))
dev.off()

idx <- (as.numeric(factor(gender_self))==1 & predicted_gender==2) | (as.numeric(factor(gender_self))==2 & predicted_gender==1)
sum(idx)
gender_self[idx]
predicted_gender[idx]

########### bad sample
De <- (RM|idx)
sum(idx)
sum(RM)
sum(De)


tab <- ID_524[De]
write.csv(tab,file="./tables/Bad_samples.csv")

########## Remove the bad sample and save RGset
RGset <- RGset[,!De]
gender_self<-gender_self[!De]
ID_slide<-ID_slide[!De]
ID_long_final<-ID_524[!De]
ID_short_524<-gender_self_order[,1]
ID_short_final<-ID_short_524[!De]
length(ID_short_final)

save(RGset,file="./RGset_QCed.rda") 
save(ID_short_final,file="./ID_short_final_QCed.rda") 

########### final preprocess
load("./RGset_QCed.rda")
object=preprocessIllumina(RGset) # preprocess
object<-mapToGenome(object)
object=ratioConvert(object, type="Illumina")#convert raw methylation data
beta <- getBeta(object) #get the beta value for each probe

dim(beta)

#### save beta
keepIndex=which(!seqnames(object)%in%c("chrX","chrY")) # mark probles on X and Y chromosome
beta <- beta[keepIndex,] #only keep probes on autosomes
dim(beta)
save(beta, file="./Beta_QCed.rda") #saves the normalized beta values to new file

## save PCs	
library(corpcor)
ss <- fast.svd(beta,tol=0) # this may take a while
save(ss,file="./fast_svd_QCed.rda")  #saves the output as a new file to be used as covariate in downstream analyses

### save Cell-type
fileName <- "./cellcount_QCed.rda"  #Output file name
cellcount <- estimateCellCounts(RGset) #estimate cell count
save(cellcount, file=fileName) #save output

dim(cellcount)
cellcount[1:3,]
################################
#######QC plots after preprocessing
################################

Cohort <- "PPMI" #Please change to your own cohort name.
mydist=dist(t(beta[sample(1:nrow(beta),10000),]))
mds=cmdscale(mydist)
#pd=pData(dat)
pdf(paste("./figures/Batch_effect_mds_",Cohort, "_QCed-4sd.pdf",sep=""),width=11,height=4)
plot(mds[,1],col=as.numeric(as.factor(ID_slide))+1,xlab="",ylab="First Principal component",xaxt="n")
dev.off()

dim(beta)
dim(mds)
length(ID_slide)
length(ID_long_final)
length(ID_short_final)
################################



