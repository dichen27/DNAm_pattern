clear;
clc;
close;

addpath(genpath('/Users/dichen/Documents/MATLAB/cd_function'));
working_dir='/Users/dichen/Documents/Methylation_Running/';
cd(working_dir);

%% load data

%t14=load('./data/t14.txt');
t_19=load('./data/t_19.mat');
t19=t_19.t19;
ID_506=t_19.t19_ID;


%% cov methylation

cov_methy_19=csvread('./data/methy_cov19_506.csv',1,0);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MDD
%% load MDD-PRS-19
% 'MDD_19_regress_ID','PRS_MDD_19_regress_ID'

load('./PRS_mental_health/(2)save_for_steiger_PRS_MDD_19.mat')


%% load  DNAm-19 cluster-1 and 7
% 'MDD_19_residulized_ID','ME_19_cluster_1_ID','ME_19_cluster_7_ID'
load('./(A-1)save_for_steiger_MDD_cluster1_7_age19.mat')


% cluster-1
ID_inter=intersect(PRS_MDD_19_regress_ID(:,1),ME_19_cluster_1_ID(:,1));

PRS_MDD_19_regress_ID_inter=cd_align_intersect(PRS_MDD_19_regress_ID,ID_inter);
ME_19_cluster_1_ID_inter=cd_align_intersect(ME_19_cluster_1_ID,ID_inter);
MDD_19_regress_ID=cd_align_intersect(MDD_19_regress_ID,ID_inter);

[r,p,df,n]=cd_partialcorr_NaN(MDD_19_regress_ID,ME_19_cluster_1_ID_inter,PRS_MDD_19_regress_ID_inter)
0.5*p


% cluster-7
ID_inter=intersect(PRS_MDD_19_regress_ID(:,1),ME_19_cluster_7_ID(:,1));

PRS_MDD_19_regress_ID_inter=cd_align_intersect(PRS_MDD_19_regress_ID,ID_inter);
ME_19_cluster_7_ID_inter=cd_align_intersect(ME_19_cluster_7_ID,ID_inter);
MDD_19_regress_ID=cd_align_intersect(MDD_19_regress_ID,ID_inter);

[r,p,df,n]=cd_partialcorr_NaN(MDD_19_regress_ID,ME_19_cluster_7_ID_inter,PRS_MDD_19_regress_ID_inter)
0.5*p

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SCZ
%% load SCZ-PRS-19
% 'SCZ_negative_19_regressed_ID','PRS_SCZ_19_regressed_ID'
load('./PRS_mental_health/(2)save_for_steiger_PRS_SCZ_19.mat')

%% load  DNAm-19 cluster-7
%'SCZ_19_residulized_ID','ME_19_cluster_7_ID')
load('./(A-1)save_for_steiger_SCZ_cluster7_age19.mat')

ID_inter=intersect(PRS_SCZ_19_regressed_ID(:,1),ME_19_cluster_7_ID(:,1));

PRS_SCZ_19_regressed_ID_inter=cd_align_intersect(PRS_SCZ_19_regressed_ID,ID_inter);
ME_19_cluster_7_ID_inter=cd_align_intersect(ME_19_cluster_7_ID,ID_inter);
SCZ_negative_19_regressed_ID_inter=cd_align_intersect(SCZ_negative_19_regressed_ID,ID_inter);


[r,p,df,n]=cd_partialcorr_NaN(SCZ_negative_19_regressed_ID_inter,ME_19_cluster_7_ID_inter,PRS_SCZ_19_regressed_ID_inter)
0.5*p
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% finished
close all; % close all figures
cd(working_dir);














