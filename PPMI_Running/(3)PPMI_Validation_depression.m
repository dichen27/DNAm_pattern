clear;
clc;
close;

addpath(genpath('/Users/dichen/Documents/LD_regression_Reply1/230211-cd_function/'));
working_dir='/Users/dichen/Documents/PPMI_Running/';
cd(working_dir);


%% load cov
load('./cov_PPMI_methy.mat')
ID_518=cov_PPMI_methy(:,1);


%% age, edu, and site
load('./cov_BL.mat')

cov_site_age_edu=[cov_BL(:,2),cov_BL(:,16),cov_BL(:,19),cov_BL(:,1)];
cov_site_age_edu=table2array(cov_site_age_edu);
cov_site_age_edu_inter=cd_align_intersect(cov_site_age_edu,ID_518);

ID_inter=cov_site_age_edu_inter(:,1);
site_inter=cd_cov_station(cov_site_age_edu_inter(:,end));
cov_1=[cov_site_age_edu_inter(:,1:end-1),site_inter]; %% age, edu, and site


%% gender and race
load('./Demographics.mat')
cov_gender_race=[Demographics(:,2),Demographics(:,10),Demographics(:,19)];
cov_gender_race=table2array(cov_gender_race);
cov_2=cd_align_intersect(cov_gender_race,ID_inter); % gender and race

cov_behavior=[cov_1,cov_2(:,2:end)];


%% cov_PPMI_methy_inter
cov_PPMI_methy_inter=cd_align_intersect(cov_PPMI_methy,ID_inter);
cov_methy=[cov_behavior,cov_PPMI_methy_inter(:,2:end)];


%% grouplabel
load('./ParticipantStatus.mat')

grouplabel=table2array(ParticipantStatus(:,1:2));
grouplabel_inter=cd_align_intersect(grouplabel,ID_inter);



%index_PD=grouplabel_inter(:,2)==1; % PD
index_HC=grouplabel_inter(:,2)==2; % HC
%index_SWEDD=grouplabel_inter(:,2)==3;% SWEDD


%% load data

%PPMI_cluster18 = readtable('./PPMI_cluster18.csv');
%PPMI_cluster18=table2array(PPMI_cluster18(:,2:end))
load('./PPMI_cluster18.mat');


PPMI_cluster18_inter=cd_align_intersect(PPMI_cluster18,ID_inter);


%% load the MDD data
load('./MDD_BL_UPDRS.mat')
MDD_BL_score_inter=cd_align_intersect(MDD_BL_score,ID_inter);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% validation in Cluster-1
PPMI_cluster_1=[PPMI_cluster18_inter(:,1:2)];

% correlation
[r,p,df,n]=cd_partialcorr_NaN(PPMI_cluster_1(index_HC,:),MDD_BL_score_inter(:,1:2),cov_methy)
0.5*p

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% validation in Cluster-7
PPMI_cluster_7=[PPMI_cluster18_inter(:,1),PPMI_cluster18_inter(:,8)];

% correlation
[r,p,df,n]=cd_partialcorr_NaN(PPMI_cluster_7(index_HC,:),MDD_BL_score_inter(:,1:2),cov_methy)
0.5*p




%% finished

close all; % close all figures
cd(working_dir);











