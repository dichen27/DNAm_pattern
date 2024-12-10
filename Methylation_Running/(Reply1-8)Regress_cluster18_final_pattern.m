clear;
clc;
close;


addpath(genpath('/Users/dichen/Documents/MATLAB/cd_function/'));
%% load data
data_dir='/Users/dichen/Documents/Methylation_Running/data/';

t14=load('/Users/dichen/Documents/Methylation_Running/data/t14.txt');
t_19=load([data_dir,'t_19.mat']);
t19=t_19.t19;
ID_506=t_19.t19_ID;


% new working_dir methylation age
working_dir='/Users/dichen/Documents/Methylation_Running/';
cd(working_dir);

%% cov 
cov_methy_14=csvread('./data/methy_cov14_506.csv',1,0);
% site wave
cov_14_con=[cov_methy_14(:,1),cov_methy_14(:,3:9),cov_methy_14(:,end-1:end)];

% cov_19
cov_methy_19=csvread('./data/methy_cov19_506.csv',1,0);
% site wave
cov_19_con=[cov_methy_19(:,1),cov_methy_19(:,3:9),cov_methy_19(:,end)];


%% load DNAm data
% age-14
data_cluster18_age14=csvread('./data_14_origin_cluster18.csv',1,1);
data_cluster18_age14_regressed=cd_regress_residual(data_cluster18_age14,cov_14_con(:,2:end));
pattern_age14_regress=corr(data_cluster18_age14_regressed,data_cluster18_age14_regressed);

save("./(Reply1-8)data_cluster18_age14_regressed.mat","data_cluster18_age14_regressed")
save("./(Reply1-8)pattern_age14_regress.mat","pattern_age14_regress")

% age-19
data_cluster18_age19=csvread('./data_19_origin_cluster18.csv',1,1);
data_cluster18_age19_regressed=cd_regress_residual(data_cluster18_age19,cov_19_con(:,2:end));
pattern_age19_regress=corr(data_cluster18_age19_regressed,data_cluster18_age19_regressed);

save("./(Reply1-8)data_cluster18_age19_regressed.mat","data_cluster18_age19_regressed")
save("./(Reply1-8)pattern_age19_regress.mat","pattern_age19_regress")

% substract
substract_cluster18_regress_19_14=data_cluster18_age19_regressed-data_cluster18_age14_regressed;

% longitudinal pattern
pattern_longitudinal_cluster18_regressed=corr(substract_cluster18_regress_19_14,substract_cluster18_regress_19_14);

save("./(Reply1-8)substract_cluster18_regress_19_14.mat","substract_cluster18_regress_19_14")
save("./(Reply1-8)pattern_longitudinal_cluster18_regressed.mat","pattern_longitudinal_cluster18_regressed")

%% finished
close all; % close all figures
cd(working_dir);




















%%











