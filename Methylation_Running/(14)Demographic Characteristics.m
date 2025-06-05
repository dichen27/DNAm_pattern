clear;
clc;
close;

addpath(genpath('/Users/dichen/Documents/LD_regression_Reply1/230211-cd_function/'));
working_dir='/Users/dichen/Documents/Methylation_Running/';
cd(working_dir);

%% load data

%t14=load('./data/t14.txt');
t_19=load('./data/t_19.mat');
t19=t_19.t19;
ID_506=t_19.t19_ID;


%% cov methylation
cov_methy_14=csvread('./data/methy_cov14_506.csv',1,0);
% cov_methy_14_inter(:,10:end)=[];

cov_methy_19=csvread('./data/methy_cov19_506.csv',1,0);
% cov_methy_19_inter(:,10:end)=[];

%% Sex
tabulate(cov_methy_14(:,2))


%% load Handness
cov_all=xlsread('./data/cov_all.xlsx');
cov_506=cd_align_intersect(cov_all,ID_506);

% Handness
tabulate(cov_506(:,end))

% Load Age
% age-14
Age_methylation_14=xlsread('./data/Age_methylation_14.xlsx');
Age_methylation_14_506=cd_align_intersect(Age_methylation_14,ID_506);

min(Age_methylation_14_506(:,end))
max(Age_methylation_14_506(:,end))

mean(Age_methylation_14_506(:,end))
std(Age_methylation_14_506(:,end))

% age-19
Age_methylation_19=xlsread('./data/Age_methylation_19.xlsx');
Age_methylation_19_506=cd_align_intersect(Age_methylation_19,ID_506);

min(Age_methylation_19_506(:,end))
max(Age_methylation_19_506(:,end))

mean(Age_methylation_19_506(:,end))
std(Age_methylation_19_506(:,end))


%% finished
close all; % close all figures
cd(working_dir);














