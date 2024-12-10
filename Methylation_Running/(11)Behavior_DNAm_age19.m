clear;
clc;
close;

addpath(genpath('/Users/dichen/Documents/MATLAB/cd_function/'));
working_dir='/Users/dichen/Documents/Methylation_Running/';
cd(working_dir);

%% load data

%t14=load('./data/t14.txt');
t_19=load('./data/t_19.mat');
t19=t_19.t19;
ID_506=t_19.t19_ID;


%% cov methylation
cov_methy_19=csvread('./data/methy_cov19_506.csv',1,0);



ME_19_all=csvread('./data_19_origin_cluster18.csv',1,1);
Module0_19=csvread('(Reply1-2)data_19_module_0.csv',1,1);
Module0_19_ID=[ID_506,Module0_19];

select_cluster=[1,3:9,11,12]';% select 10 clusters
length(select_cluster)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Substance Use with DNAm clusters
data_use_dir='./IMAGEN_substance_total_FU2.xlsx';
data_all=importdata(data_use_dir);
phenotype_19=data_all.data;



ME_19_select_ID=[ID_506,ME_19_all(:,select_cluster)];


[r,P_all,df,n]=cd_partialcorr_NaN(phenotype_19,ME_19_select_ID,cov_methy_19);
P_fdr=cd_FDR(P_all);
data_all.colheaders


%% intersect for controlling
ID_inter=intersect(ME_19_select_ID(:,1),phenotype_19(:,1));
phenotype_19_inter=cd_align_intersect(phenotype_19,ID_inter);
cov_methy_19_inter=cd_align_intersect(cov_methy_19,ID_inter);
ME_19_common_ID=[ME_19_select_ID(:,1),ME_19_select_ID(:,6),ME_19_select_ID(:,9:11)];

%% DNAm corr smoking: controlling Marijuana
smoking_19=[phenotype_19_inter(:,1),phenotype_19_inter(:,3)];
cov_final=[cov_methy_19_inter,phenotype_19_inter(:,5)];% add Marijuana as cov

[r,p,df,n]=cd_partialcorr_NaN(smoking_19,ME_19_common_ID,cov_final)



%% DNAm corr Marijuana: controlling smoking
Marijuana_19=[phenotype_19_inter(:,1),phenotype_19_inter(:,5)];
cov_final=[cov_methy_19_inter,phenotype_19_inter(:,3)];% add smoking as cov

[r,p,df,n]=cd_partialcorr_NaN(Marijuana_19,ME_19_common_ID,cov_final)



%% subtance use module-0
data_use_dir='./IMAGEN_substance_total_FU2.xlsx';
data_all=importdata(data_use_dir);
phenotype_19=data_all.data;

[r,p_subtance,df,n]=cd_partialcorr_NaN(Module0_19_ID,phenotype_19,cov_methy_19)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MDD module-0
data_use_dir='./FU2_Behaviour/IMAGEN_FU2_ADRS.csv';
data_all=importdata(data_use_dir);
phenotype_19=data_all.data(:,1:2);

[r,p_MDD,df,n]=cd_partialcorr_NaN(Module0_19_ID,phenotype_19,cov_methy_19)
P_module0_19=[p_subtance';p_MDD];


%% MDD with DNAm clusters
data_use_dir='./FU2_Behaviour/IMAGEN_FU2_ADRS.csv';
data_all=importdata(data_use_dir);
phenotype_19=data_all.data(:,1:2);
ME_19_select_ID=[ID_506,ME_19_all(:,select_cluster)];

% main result: corr
[r,P_all,df,n]=cd_partialcorr_NaN(phenotype_19,ME_19_select_ID,cov_methy_19)

% FDR-correction
P_dfr=cd_FDR(P_all)

%% save for steiger
ID_inter=intersect(phenotype_19(:,1),ME_19_select_ID(:,1));
ID_inter=intersect(ID_inter,cov_methy_19(:,1));

cov_methy_19_inter=cd_align_intersect(cov_methy_19,ID_inter);
phenotype_19_inter=cd_align_intersect(phenotype_19,ID_inter);
ME_19_select_ID=cd_align_intersect(ME_19_select_ID,ID_inter);

% residual for steiger
MDD_19_residulized=cd_regress_residual(phenotype_19_inter(:,2),cov_methy_19_inter(:,2:9));% sex,site

ME_19_cluster_1=cd_regress_residual(ME_19_select_ID(:,2),cov_methy_19_inter(:,2:end));
ME_19_cluster_7=cd_regress_residual(ME_19_select_ID(:,7),cov_methy_19_inter(:,2:end));


MDD_19_residulized_ID=[ID_inter,MDD_19_residulized];
ME_19_cluster_1_ID=[ID_inter,ME_19_cluster_1];
ME_19_cluster_7_ID=[ID_inter,ME_19_cluster_7];

% save for steiger
save('./(A-1)save_for_steiger_MDD_cluster1_7_age19.mat','MDD_19_residulized_ID','ME_19_cluster_1_ID','ME_19_cluster_7_ID')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SCZ with DNAm clusters
data_use_dir='./FU2_Behaviour/IMAGEN_FU2_CAPE.csv';
data_all=importdata(data_use_dir);
phenotype_19=data_all.data;
data_all.textdata

ME_19_select_ID=[ID_506,ME_19_all(:,select_cluster)];


[r,P_all,df,n]=cd_partialcorr_NaN(phenotype_19,ME_19_select_ID,cov_methy_19)
P_postive=P_all(3,:);
P_negetive=P_all(9,:);
P_depressive=P_all(6,:);

ID_inter=intersect(phenotype_19(:,1),ME_19_select_ID(:,1));
ID_inter=intersect(ID_inter,cov_methy_19(:,1));

cov_methy_19_inter=cd_align_intersect(cov_methy_19,ID_inter);
phenotype_19_inter=cd_align_intersect(phenotype_19,ID_inter);
ME_19_select_ID=cd_align_intersect(ME_19_select_ID,ID_inter);

% residual
SCZ_19_residulized=cd_regress_residual(phenotype_19_inter(:,10),cov_methy_19_inter(:,2:9));
ME_19_cluster_7=cd_regress_residual(ME_19_select_ID(:,7),cov_methy_19_inter(:,2:end));


SCZ_19_residulized_ID=[ID_inter,SCZ_19_residulized];
ME_19_cluster_7_ID=[ID_inter,ME_19_cluster_7];

%% save for steiger
save('./(A-1)save_for_steiger_SCZ_cluster7_age19.mat','SCZ_19_residulized_ID','ME_19_cluster_7_ID')


%% finished
close all; % close all figures
cd(working_dir);














