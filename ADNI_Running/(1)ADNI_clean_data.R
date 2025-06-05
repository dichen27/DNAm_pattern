rm(list=ls())
library(R.matlab) 
library(xlsx)
setwd("/Users/dichen/Documents/ADNI_Running/")
##########################
cd_maxmin_date <- function(date_all){
  date_order<-matrix(data=NA, nrow = length(date_all), ncol =2, byrow = FALSE, dimnames = NULL)
  for (i in 1:length(date_all)) {
    
    index_max<-date_all[i]>=date_all
    date_order[i,1]<-sum(index_max)
    
    index_min<-date_all[i]<=date_all
    date_order[i,2]<-sum(index_min)
  }
  
  index_max<-date_order[,1]==length(date_all)
  date_max<-date_all[index_max]
  
  index_min<-date_order[,2]==length(date_all)
  date_min<-date_all[index_min]
  
  return (c(date_max[1],date_min[1]))
}

##########################
## load the ADNI data
data_dir="/Users/dichen/Documents/ADNI_Running/ADNI/data_QCed/"
load(paste(data_dir,"Beta_QCed.rda",sep=""))
row.names(beta) [1]
length(row.names(beta) )
cpg_name_ADNI<-row.names(beta)


data_ADNI<-t(beta)
row.names(beta) [1]

dim(data_ADNI)
ID_long_final<-row.names(data_ADNI) 
ID_long_final[3]

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
cov_ADNI_methy<-cbind(ID_long_final,ID_short_final,cellcount,PCs)
write.csv(cov_ADNI_methy, file = "cov_ADNI_methy.csv")

dim(cov_ADNI_methy)
cov_ADNI_methy[1:3,]
################ Grouping
Methylation_SampleAnnotation<-read.csv('./ADNI_DNA_Methylation_SampleAnnotation_20170530.csv')

Methylation_SampleAnnotation$Edate = as.Date(Methylation_SampleAnnotation$Edate)
Methylation_SampleAnnotation$DateDrawn=as.Date(Methylation_SampleAnnotation$DateDrawn)

# align the 1889
Methylation_SampleAnnotation_1889<-data.frame()
for (i in 1:length(ID_long_final)) {

  index<-Methylation_SampleAnnotation$barcodes==ID_long_final[i]

  Methylation_SampleAnnotation_1889<-rbind( Methylation_SampleAnnotation_1889,Methylation_SampleAnnotation[index,])
}


# grouping
index_ADNI_GO<-Methylation_SampleAnnotation_1889[,3]=="ADNIGO";
index_ADNI_2<-Methylation_SampleAnnotation_1889[,3]=="ADNI2";

sum(index_ADNI_GO)
sum(index_ADNI_2)

ID_long_ADNI_GO<-Methylation_SampleAnnotation_1889$barcodes[index_ADNI_GO]
ID_long_ADNI_2<-Methylation_SampleAnnotation_1889$barcodes[index_ADNI_2]


ID_short_ADNI_GO<-Methylation_SampleAnnotation_1889$RID[index_ADNI_GO]
ID_short_ADNI_2<-Methylation_SampleAnnotation_1889$RID[index_ADNI_2]


length(ID_short_ADNI_GO)
length(ID_short_ADNI_2)
length(unique(ID_short_ADNI_2))


##################### ID_long_ADNI_2_unique
ID_short_ADNI_2_unique<-unique(ID_short_ADNI_2)
ID_long_ADNI_2_unique<-matrix(data=NA, nrow = length(ID_short_ADNI_2_unique), ncol =1, byrow = FALSE, dimnames = NULL)
date_max_ADNI_2_unique<-Methylation_SampleAnnotation_1889$DateDrawn[1]
for (i in 1:length(ID_short_ADNI_2_unique)) {
#i<-277
index<-Methylation_SampleAnnotation_1889$RID==ID_short_ADNI_2_unique[i]
ID_long_i<-Methylation_SampleAnnotation_1889$barcodes[index]

if (sum(is.na(Methylation_SampleAnnotation_1889$Edate[index]))==sum(index)){ID_long_ADNI_2_unique[i]<-ID_long_i}
else{
date_all<-Methylation_SampleAnnotation_1889$Edate[index]
date_max<-cd_maxmin_date(date_all)[1]

index_max<-date_all==date_max
ID_long_max<-ID_long_i[index_max]

date_max_ADNI_2_unique[i]<-date_max
ID_long_ADNI_2_unique[i]<-ID_long_max[1]
}
}

length(ID_long_ADNI_2_unique)

##########################################################
##################### Mapping
grouplabel_information<-read.csv('./PTDEMOG.csv')


grouplabel_information_mapped<-data.frame()
for (i in 1:length(ID_short_ADNI_2_unique)) {
  i
 
  index_1<-(grouplabel_information$RID==ID_short_ADNI_2_unique[i] & grouplabel_information$Phase=="ADNI2")# RID, "ADNI2"
  sum(index_1)
  if (sum(index_1)>0){
    grouplabel_information_select_1<-grouplabel_information[index_1,]
    
    # dele the row with all NA
    index<-is.na(grouplabel_information_select_1$PTCOGBEG)&is.na(grouplabel_information_select_1$PTADDX)
    grouplabel_information_select<-grouplabel_information_select_1[!index, ] 
    
    
    num_AD<-sum(grouplabel_information_select$PTADDX<9900 & !is.na(grouplabel_information_select$PTADDX))
    num_MCI<-sum((grouplabel_information_select$PTCOGBEG<9900&!is.na(grouplabel_information_select$PTCOGBEG)) & (grouplabel_information_select$PTADDX==9999|is.na(grouplabel_information_select$PTADDX)))
    num_HC<-sum((grouplabel_information_select$PTCOGBEG==9999|is.na(grouplabel_information_select$PTCOGBEG)) & grouplabel_information_select$PTADDX==9999)
  
    #num_AD<-sum(grouplabel_information_select$PTCOGBEG<9900 & grouplabel_information_select$PTADDX<9900 & !is.na(grouplabel_information_select$PTADDX))
    #num_MCI<-sum(grouplabel_information_select$PTCOGBEG<9900 & (grouplabel_information_select$PTADDX==9999|is.na(grouplabel_information_select$PTADDX)))
    #num_HC<-sum((grouplabel_information_select$PTCOGBEG==9999|is.na(grouplabel_information_select$PTCOGBEG)) & grouplabel_information_select$PTADDX==9999)
    
    if (num_AD==1){grouplabel=2}
    if (num_MCI==1){grouplabel=1}
    if (num_HC==1){grouplabel=0}
    
    
    RID<-grouplabel_information_select$RID
    viscode<-grouplabel_information_select$VISCODE
    viscode2<-grouplabel_information_select$VISCODE2
    
    site<-grouplabel_information_select$SITEID
    gender<-grouplabel_information_select$PTGENDER
    educ<-grouplabel_information_select$PTEDUCAT
    
    grouplabel_information_mapped<-rbind(grouplabel_information_mapped,c(ID_long_ADNI_2_unique[i],RID,viscode,viscode2,site,gender,educ,grouplabel))
  }
}
names(grouplabel_information_mapped)<-c("ID_long_ADNI_2_unique","RID","viscode","viscode2","site","gender","educ","grouplabel")


dim(grouplabel_information_mapped)
grouplabel_information_mapped[1:3,]

##################################### age
age_information<-read.csv('./ADNI_age.csv')
age_all<-age_information$AGE

age_mapped<-data.frame()
for (i in 1:length(ID_short_ADNI_2_unique)) {
 # i<-122
  index<-age_information$RID==ID_short_ADNI_2_unique[i] 
  sum(index)
 
  age_all<-age_information$AGE[index]
  age_all
  
  
  # Use the all() function to check if all elements in the vector are equal
  if (all(age_all == age_all[1])) {
    age_mapped<-rbind(age_mapped,c(ID_short_ADNI_2_unique[i],age_all[1]))
  } else {
    cat("Not all elements are equal\n")
  }
}
names(age_mapped)<-c("ID_short_ADNI_2_unique","age")
mean(age_mapped$age)
sd(age_mapped$age)

writeMat("./age_606.mat",age_606=age_mapped)
########################################################## align and save 606 
# save data
data_ADNI_606<-matrix(data=NA, nrow = 606, ncol =dim(data_ADNI)[2], byrow = FALSE, dimnames = NULL)
for (i in 1:606) {
print(i)
  
index<-ID_long_final==ID_long_ADNI_2_unique[i]
data_ADNI_606[i,]<-data_ADNI[index,]
}
save(data_ADNI_606,file="./data_ADNI2_606/data_ADNI_606.rda") 


# save cov
cov_ADNI_methy_606<-matrix(data=NA, nrow = 606, ncol =dim(cov_ADNI_methy)[2], byrow = FALSE, dimnames = NULL)
for (i in 1:606) {
print(i)
  
index<-cov_ADNI_methy[,1]==ID_long_ADNI_2_unique[i]
cov_ADNI_methy_606[i,]<-cov_ADNI_methy[index,]
}
save(cov_ADNI_methy_606, file = "./data_ADNI2_606/cov_ADNI_methy_606.rda")


#save ID
save(ID_long_ADNI_2_unique, file = "./data_ADNI2_606/ID_long_ADNI_606.rda")
save(ID_short_ADNI_2_unique, file = "./data_ADNI2_606/ID_short_ADNI_606.rda")
save(cpg_name_ADNI, file = "./data_ADNI2_606/cpg_name_ADNI.rda")


grouplabel_606<-grouplabel_information_mapped
save(grouplabel_606, file = "./data_ADNI2_606/grouplabel_606.rda")

grouplabel_606[1:3,]
tabulate(factor(grouplabel_606$gender))
tabulate(factor(grouplabel_606$educ))
tabulate(factor(grouplabel_606$grouplabel))
##########################################################
#####################finished








