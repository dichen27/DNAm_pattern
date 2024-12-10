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
cov_methy_14=csvread('./data/methy_cov14_506.csv',1,0);




%% load PRS 

load('./PRS_mental_health/(2)save_for_steiger_PRS_MDD_14.mat')
% 'MDD_14_regress_ID','PRS_MDD_14_regress_ID')
% 'ME_14_residulized_cluster_7_ID','MDD_14_residulized'



%% load  DNAm cluster-7
load('./(A-2)save_for_steiger_MDD_cluster7_14.mat')

[r,p,df,n]=cd_partialcorr_NaN(MDD_14_regress_ID,ME_14_residulized_cluster_7_ID,PRS_MDD_14_regress_ID)
0.5*p


%% finished
close all; % close all figures
cd(working_dir);














