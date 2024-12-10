clear;
clc;
close;

addpath(genpath('/Users/dichen/Documents/MATLAB/cd_function/'));
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

tabulate(grouplabel_inter(:,2))


load('./PPMI_ME175.mat');
PPMI_ME175_inter=cd_align_intersect(PPMI_ME175,ID_inter);


%% grouping

cov_site=[cov_1(:,1),cov_1(:,4:end)];
cov_site_inter=cd_align_intersect(cov_site,ID_inter);


Methy_ME175_site_regressed=cd_regress_residual(PPMI_ME175_inter(:,2:end),cov_site_inter(:,2:end));



index_AD=grouplabel_inter(:,2)==1;
sum(index_AD)
ID_AD=ID_inter(index_AD);



index_HC=grouplabel_inter(:,2)==2;
sum(index_HC)
ID_HC=ID_inter(index_HC);



index_SWEDD=grouplabel_inter(:,2)==3;
sum(index_SWEDD)
ID_SWEDD=ID_inter(index_SWEDD);


save('./(2-2)ID_grouping.mat','ID_AD','ID_HC','ID_SWEDD')



%% finished
close all; % close all figures
cd(working_dir);











