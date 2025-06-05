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



index_1=cov_14_con(:,end-1)==1;
sum(index_1)

index_2=cov_14_con(:,end)==1;
sum(index_2)

index_850k=cov_14_con(:,end)==0&cov_14_con(:,end-1)==0;
sum(index_850k)

index_450k=cov_14_con(:,end)==1|cov_14_con(:,end-1)==1;
sum(index_450k)

t14_450k=t14(index_450k,:);
t14_850k=t14(index_850k,:);

save([working_dir,'index_450k_850k.mat'],'index_450k','index_850k');




cov_14_450k=cov_14_con(index_450k,:);
cov_14_450k(:,end)=[];

cov_14_850k=cov_14_con(index_850k,:);
cov_14_850k(:,end-1:end)=[];


% cov_19
cov_methy_19=csvread('./data/methy_cov19_506.csv',1,0);

% site wave
cov_19_con=[cov_methy_19(:,1),cov_methy_19(:,3:9),cov_methy_19(:,end)];


%% regress respectively

% age-14
t14_450k_regress=cd_regress_residual(t14_450k,cov_14_450k(:,2:end));
t14_450k_regress=[ID_506(index_450k),t14_450k_regress];
save([working_dir,'data_regressed/t14_450k_regress.mat'],'t14_450k_regress');


t14_850k_regress=cd_regress_residual(t14_850k,cov_14_850k(:,2:end));
t14_850k_regress=[ID_506(index_850k),t14_850k_regress];
save([working_dir,'data_regressed/t14_850k_regress.mat'],'t14_850k_regress');


% age-19
t19_regress=cd_regress_residual(t19,cov_19_con(:,2:end));


t19_450k_regress=t19_regress(index_450k,:);
t19_450k_regress=[ID_506(index_450k),t19_450k_regress];
save([working_dir,'data_regressed/t19_450k_regress.mat'],'t19_450k_regress');


t19_850k_regress=t19_regress(index_850k,:);
t19_850k_regress=[ID_506(index_850k),t19_850k_regress];
save([working_dir,'data_regressed/t19_850k_regress.mat'],'t19_850k_regress');
%% regress all 506 individual

% age-14
t14_full_regress=cd_regress_residual(t14,cov_14_con(:,2:end));
t14_full_regress=[ID_506,t14_full_regress];
save([working_dir,'data_regressed/t14_full_regress.mat'],'t14_full_regress');

% age-19
t19_full_regress=cd_regress_residual(t19,cov_19_con(:,2:end));
t19_full_regress=[ID_506,t19_full_regress];
save([working_dir,'data_regressed/t19_full_regress.mat'],'t19_full_regress');


%% finished
close all; % close all figures
cd(working_dir);




















%%











