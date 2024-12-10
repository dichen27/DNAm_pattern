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

%% load T1
load('./(9)Freesurfer_data.mat');
%'GrayVol_14','GrayVol_19','SurfArea_14','SurfArea_19','ThickAvg_14','ThickAvg_19','SubcortexVol_14','SubcortexVol_19','ICV_14','ICV_19'


%%%%%%%%%%% 
%% The following four lines should be muted in sequence, and then saved separately.
ID_inter=intersect(SurfArea_14(:,1),SurfArea_19(:,1));
 % ID_inter=intersect(ThickAvg_14(:,1),ThickAvg_19(:,1));
  % ID_inter=intersect(GrayVol_14(:,1),GrayVol_19(:,1));
 %ID_inter=intersect(SubcortexVol_14(:,1),SubcortexVol_19(:,1));


ID_inter=intersect(ID_inter,ID_506);

T1_14_inter=cd_align_intersect(SubcortexVol_14,ID_inter);
T1_19_inter=cd_align_intersect(SubcortexVol_19,ID_inter);

T1_change_inter=T1_19_inter(:,2:end)-T1_14_inter(:,2:end);
%T1_change_inter(:,33:end)=[];
T1_change_ID=[ID_inter,T1_change_inter];

%% load TIV
ICV_14_inter=cd_align_intersect(ICV_14,ID_inter);
ICV_19_inter=cd_align_intersect(ICV_19,ID_inter);



%% load ME
%ME_second_14_all=csvread('./ME_second_Cluster18_weight_BL.csv',1,1);
ME_second_14_all=csvread('./data_14_origin_cluster18.csv',1,1);


ME_second_14_all_ID=[ID_506,ME_second_14_all];
ME_14_inter=cd_align_intersect(ME_second_14_all_ID,ID_inter);


%ME_second_19_all=csvread('./ME_second_Cluster18_weight_FU2.csv',1,1);
ME_second_19_all=csvread('./data_19_origin_cluster18.csv',1,1);


ME_second_19_all_ID=[ID_506,ME_second_19_all];
ME_19_inter=cd_align_intersect(ME_second_19_all_ID,ID_inter);

ME_change=ME_19_inter(:,2:end)-ME_14_inter(:,2:end);
mean(ME_change)

% [h,p,ci,stats] = ttest2(ME_19_inter(:,2),ME_14_inter(:,2))
ME_change_ID=[ID_inter,ME_change];

%% cov_no_methylation_468
SDQ14=xlsread('./data/sdq14_ori_cov.xlsx');
cov_gender_site=[SDQ14(:,1),SDQ14(:,4:end)];

cov_no_methylation_468=cd_align_intersect(cov_gender_site,ID_inter);

cov_no_methylation_TIV_14=[cov_no_methylation_468,ICV_14_inter(:,2)];
cov_no_methylation_TIV_19=[cov_no_methylation_468,ICV_19_inter(:,2)];

%% regressed ME
cov_methy_14=csvread('./data/methy_cov14_506.csv',1,0);
cov_methy_14_inter=cd_align_intersect(cov_methy_14,ID_inter);

wave_14=cov_methy_14_inter(:,end-1:end);

index_1=cov_methy_14_inter(:,end-1)==1;
sum(index_1)

index_2=cov_methy_14_inter(:,end)==1;
sum(index_2)

index_3=cov_methy_14_inter(:,end)==0&cov_methy_14_inter(:,end-1)==0;
sum(index_3)

index_4=cov_methy_14_inter(:,end)==1|cov_methy_14_inter(:,end-1)==1;
sum(index_4)

% cov_methy_14_inter(:,10:end)=[];

cov_methy_19=csvread('./data/methy_cov19_506.csv',1,0);
cov_methy_19_inter=cd_align_intersect(cov_methy_19,ID_inter);
wave_19=cov_methy_19_inter(:,end);
% cov_methy_19_inter(:,10:end)=[];


%% cov_final
cov_14=[cov_methy_14_inter,ICV_14_inter(:,2)];
cov_19=[cov_methy_19_inter,ICV_19_inter(:,2)];

%% regressed ME
ME_14_regressed=cd_regress_residual(ME_14_inter(:,2:end),cov_methy_14_inter(:,2:end));
ME_19_regressed=cd_regress_residual(ME_19_inter(:,2:end),cov_methy_19_inter(:,2:end));

ME_14_regressed=cd_regress_residual(ME_14_inter(:,2:end),cov_14(:,2:end));
ME_19_regressed=cd_regress_residual(ME_19_inter(:,2:end),cov_19(:,2:end));


ME_change_regressed=ME_19_regressed-ME_14_regressed;


%% regressed T1

T1_14_regressed=cd_regress_residual(T1_14_inter(:,2:end),cov_no_methylation_TIV_14(:,2:end));
T1_19_regressed=cd_regress_residual(T1_19_inter(:,2:end),cov_no_methylation_TIV_19(:,2:end));

T1_14_regressed=cd_regress_residual(T1_14_inter(:,2:end),cov_14(:,2:end));
T1_19_regressed=cd_regress_residual(T1_19_inter(:,2:end),cov_19(:,2:end));

T1_change_regressed=T1_19_regressed-T1_14_regressed;



%% save data

save('./data_CCA_change_Methy_SurfArea.mat','ME_change_regressed','T1_change_regressed')
% save('./data_CCA_change_Methy_ThickAvg.mat','ME_change_regressed','T1_change_regressed')
%save('./data_CCA_change_Methy_GrayVol.mat','ME_change_regressed','T1_change_regressed')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% re-order subcortex
ID_ROI_16=[15;13;6;12;8;7;5;1;29;28;24;27;26;25;23;19];% order by plot enigma

T1_change_regressed_16=[];
for i=1:length(ID_ROI_16)
n=ID_ROI_16(i)

lin=T1_change_regressed(:,n);

T1_change_regressed_16=[T1_change_regressed_16,lin];

end

T1_change_regressed=T1_change_regressed_16;

save('./data_CCA_change_Methy_SubcortexVol.mat','ME_change_regressed','T1_change_regressed')




ID_467=ID_inter;
save('./ID_467.mat','ID_467')

%% finished

close all; % close all figures
cd(working_dir);







