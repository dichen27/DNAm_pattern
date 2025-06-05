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


%% load ME_175

ME_14_all=csvread('./data_14_origin_ME175.csv',1,1);


%% plot MDD-band-14
ME_14_all_ID=[ID_506,ME_14_all];

data_all=importdata('./IMAGEN_dawba_cleaned_BL.xlsx');
phenotype_14=data_all.data(:,1:2);

[r,p,df,n]=cd_partialcorr_NaN(ME_14_all_ID,phenotype_14,cov_methy_14)
output_r_p=[r,p]


%% reorder by cluster
% load matrix_Module_Cluster
load("(5)matrix_Module_Cluster_percent_10.mat")

matrix_Module_Cluster=double(matrix_Module_Cluster);
output=[matrix_Module_Cluster,output_r_p];

output_sort = sortrows(output, 2);



%% finished
close all; % close all figures
cd(working_dir);














