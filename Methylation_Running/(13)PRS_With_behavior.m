clear;
clc;
close;

addpath(genpath('/Users/dichen/Documents/LD_regression_Reply1/230211-cd_function/'));
working_dir='/Users/dichen/Documents/Methylation_Running/PRS_mental_health/';
cd(working_dir);

%% load data

%t14=load('./data/t14.txt');
t_19=load('../data/t_19.mat');
t19=t_19.t19;
ID_506=t_19.t19_ID;


%% cov methylation
cov_methy_14=csvread('../data/methy_cov14_506.csv',1,0);
cov_methy_19=csvread('../data/methy_cov19_506.csv',1,0);




%% PRS mental health
PRS=load('./PRS_IMAGEN_mental_health.mat');
PRS_mental_health=PRS.PRS_IMAGEN_mental_health;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MDD-14 with PRS-MDD
PRS_MDD=[PRS_mental_health(:,1),PRS_mental_health(:,2)];

data_all=importdata('/Users/dichen/Documents/Methylation_Running/PRS_mental_health/data/IMAGEN_dawba_MDD_BL.xlsx');
phenotype=data_all.data;


[r,p,df,n]=cd_partialcorr_NaN(phenotype,PRS_MDD,cov_methy_14(:,1:9))
0.5*p

%% save for steiger
ID_intersect=intersect(phenotype(:,1),cov_methy_14(:,1));
ID_intersect=intersect(PRS_MDD(:,1),ID_intersect);

phenotype_inter=cd_align_intersect(phenotype,ID_intersect);
PRS_MDD_inter=cd_align_intersect(PRS_MDD,ID_intersect);
cov_methy_14_inter=cd_align_intersect(cov_methy_14,ID_intersect);

% find n
length(cd_dele_nan_row(phenotype_inter))

% regress
MDD_14_regress=cd_regress_residual(phenotype_inter(:,2),cov_methy_14_inter(:,2:9));
PRS_MDD_14_regress=cd_regress_residual(PRS_MDD_inter(:,2),cov_methy_14_inter(:,2:9));


MDD_14_regress_ID=[ID_intersect,MDD_14_regress];
PRS_MDD_14_regress_ID=[ID_intersect,PRS_MDD_14_regress];

% save for steiger
save('./(2)save_for_steiger_PRS_MDD_14.mat','MDD_14_regress_ID','PRS_MDD_14_regress_ID')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MDD-19 with PRS-MDD
PRS_MDD=[PRS_mental_health(:,1),PRS_mental_health(:,2)];


data_all=importdata('/Users/dichen/Documents/Methylation_Running/PRS_mental_health/data/IMAGEN_FU2_ADRS.csv');
phenotype=data_all.data;


[r,p,df,n]=cd_partialcorr_NaN(phenotype,PRS_MDD,cov_methy_19(:,1:9))
0.5*p

%% save for steiger
ID_intersect=intersect(phenotype(:,1),cov_methy_19(:,1));
ID_intersect=intersect(PRS_MDD(:,1),ID_intersect);

phenotype_inter=cd_align_intersect(phenotype,ID_intersect);
PRS_MDD_inter=cd_align_intersect(PRS_MDD,ID_intersect);
cov_methy_19_inter=cd_align_intersect(cov_methy_19,ID_intersect);

% find n
length(cd_dele_nan_row(phenotype_inter))

% regress
MDD_19_regress=cd_regress_residual(phenotype_inter(:,2),cov_methy_19_inter(:,2:9));
PRS_MDD_19_regress=cd_regress_residual(PRS_MDD_inter(:,2),cov_methy_19_inter(:,2:9));

MDD_19_regress_ID=[ID_intersect,MDD_19_regress];
PRS_MDD_19_regress_ID=[ID_intersect,PRS_MDD_19_regress];

% save for steiger
save('./(2)save_for_steiger_PRS_MDD_19.mat','MDD_19_regress_ID','PRS_MDD_19_regress_ID')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SCZ-positive with PRS-SCZ

PRS_SCZ=[PRS_mental_health(:,1),PRS_mental_health(:,7)];
data_all=importdata('/Users/dichen/Documents/Methylation_Running/PRS_mental_health/data/IMAGEN_positive_SCZ_FU2.xlsx');
phenotype=data_all.data;

[r,p,df,n]=cd_partialcorr_NaN(phenotype,PRS_SCZ,cov_methy_14(:,1:9))
P_SCZ=0.5*p

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SCZ-negetive with PRS-SCZ
PRS_SCZ=[PRS_mental_health(:,1),PRS_mental_health(:,7)];
data_all=importdata('/Users/dichen/Documents/Methylation_Running/PRS_mental_health/data/IMAGEN_negetive_SCZ_FU2.xlsx');
phenotype=data_all.data;


[r,p,df,n]=cd_partialcorr_NaN(phenotype,PRS_SCZ,cov_methy_14(:,1:9))
P_SCZ=[P_SCZ;0.5*p]

%% save for steiger
ID_intersect=intersect(phenotype(:,1),cov_methy_19(:,1));
ID_intersect=intersect(PRS_SCZ(:,1),ID_intersect);

phenotype_inter=cd_align_intersect(phenotype,ID_intersect);
PRS_SCZ_inter=cd_align_intersect(PRS_SCZ,ID_intersect);
cov_methy_19_inter=cd_align_intersect(cov_methy_19,ID_intersect);

%% regress
SCZ_negative_19_regressed=cd_regress_residual(phenotype_inter(:,2),cov_methy_19_inter(:,2:9));
PRS_SCZ_19_regressed=cd_regress_residual(PRS_SCZ_inter(:,2),cov_methy_19_inter(:,2:9));

SCZ_negative_19_regressed_ID=[ID_intersect,SCZ_negative_19_regressed];
PRS_SCZ_19_regressed_ID=[ID_intersect,PRS_SCZ_19_regressed];

% save for steiger
save('./(2)save_for_steiger_PRS_SCZ_19.mat','SCZ_negative_19_regressed_ID','PRS_SCZ_19_regressed_ID')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SCZ-depress with PRS-SCZ

PRS_SCZ=[PRS_mental_health(:,1),PRS_mental_health(:,7)];
data_all=importdata('/Users/dichen/Documents/Methylation_Running/IMAGEN_FU2_SCZ_depressive_total.csv');
phenotype=data_all.data;


[r,p,df,n]=cd_partialcorr_NaN(phenotype,PRS_SCZ,cov_methy_14(:,1:9))
P_SCZ=[P_SCZ;0.5*p];


%% FDR correction
P_SCZ_FDR=cd_FDR(P_SCZ)

%% finished
close all; % close all figures
cd(working_dir);














