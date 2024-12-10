clear;
clc;
close;

addpath(genpath('/Users/dichen/Documents/MATLAB/cd_function/'));
working_dir='/Users/dichen/Documents/ADNI_Running/';
cd(working_dir);


%% load cov


%cov_ADNI_methy_606=covADNImethy606;
%cov_ADNI_methy_606(1,:)=[];
%cov_ADNI_methy_606=table2array(cov_ADNI_methy_606);
%save('./cov_ADNI_methy_606.mat','cov_ADNI_methy_606')

load('./cov_ADNI_methy_606.mat');


%% age, edu, and site
%grouplabel_606=grouplabel606
%grouplabel_606(1,:)=[];
%grouplabel_606=table2array(grouplabel_606);% site gender educ grouplabel:(2=AD 1=MCI 0 =HC)
%save('./grouplabel_606.mat','grouplabel_606')




load('./grouplabel_606.mat');% ID site gender educ grouplabel:(2=AD 1=MCI 0 =HC)

site_source=grouplabel_606(:,2);
site_matrix=cd_cov_station(site_source);

%gender
tabulate(grouplabel_606(:,3))


%age

%% load ME-175 data

%ADNIME175(:,1:2)=[];
%ADNI_ME175=table2array(ADNIME175);
%save('./ADNI_ME175.mat','ADNI_ME175');

load('./ADNI_ME175.mat')



%regrees
ADNI_ME175_site_regressed=cd_regress_residual(ADNI_ME175(:,2:end),site_matrix);

%% grouping

grouplabel=grouplabel_606(:,end);% 2=AD 1=MCI 0 =HC
tabulate(grouplabel)



index_AD=grouplabel==2;
sum(index_AD)
ID_AD=ADNI_ME175(index_AD,1);
Methy_ME175_regressed_AD=ADNI_ME175_site_regressed(index_AD,:);


index_HC=grouplabel==0;
sum(index_HC)
ID_HC=ADNI_ME175(index_HC,1);
Methy_ME175_regressed_HC=ADNI_ME175_site_regressed(index_HC,:);


index_MCI=grouplabel==1;
sum(index_MCI)
ID_MCI=ADNI_ME175(index_MCI,1);
Methy_ME175_regressed_MCI=ADNI_ME175_site_regressed(index_MCI,:);


%% AD
% ID site gender educ grouplabel:(2=AD 1=MCI 0 =HC)
cov_AD=cd_align_intersect(grouplabel_606,ID_AD);

% AD gender
tabulate(cov_AD(:,3))

% AD education
mean(cov_AD(:,4))
std(cov_AD(:,4))



%% MCI
cov_MCI=cd_align_intersect(grouplabel_606,ID_MCI);

% AD gender
tabulate(cov_MCI(:,3))

% AD education
mean(cov_MCI(:,4))
std(cov_MCI(:,4))



%% HC
% ID site gender educ grouplabel:(2=AD 1=MCI 0 =HC)
cov_HC=cd_align_intersect(grouplabel_606,ID_HC);

% HC gender
tabulate(cov_HC(:,3))

% HC education
mean(cov_HC(:,4))
std(cov_HC(:,4))



%% Age
load('./age_606.mat')
age_606_all=[age_606.ID_short_ADNI_2_unique,age_606.age];

% age AD
age_AD=cd_align_intersect(age_606_all,ID_AD);
mean(age_AD(:,2))
std(age_AD(:,2))


% age MCI
age_MCI=cd_align_intersect(age_606_all,ID_MCI);
mean(age_MCI(:,2))
std(age_MCI(:,2))


% age HC
age_HC=cd_align_intersect(age_606_all,ID_HC);
mean(age_HC(:,2))
std(age_HC(:,2))



%% finished
close all; % close all figures
cd(working_dir);











