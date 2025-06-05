clear;
clc;
close;

addpath(genpath('/Users/dichen/Documents/MATLAB/cd_function/'));
working_dir='/Users/dichen/Documents/ADNI_Running/';
cd(working_dir);

load('./cov_ADNI_methy_606.mat');


load('./grouplabel_606.mat');% ID site gender educ grouplabel:(2=AD 1=MCI 0 =HC)

site_source=grouplabel_606(:,2);
site_matrix=cd_cov_station(site_source);

%gender
tabulate(grouplabel_606(:,3))


%age

%% load ME-175 data
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
save('./(3)ADNI_ID_grouping.mat','ID_AD','ID_HC','ID_MCI')

%save('./(3)ADNI_Methy_ME175_regressed_grouping.mat','Methy_ME175_regressed_AD','Methy_ME175_regressed_HC','Methy_ME175_regressed_MCI')






%% finished
close all; % close all figures
cd(working_dir);











