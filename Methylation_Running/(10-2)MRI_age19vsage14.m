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
%% The following four lines should be muted in sequence separately.
 %ID_inter=intersect(SurfArea_14(:,1),SurfArea_19(:,1));
 % ID_inter=intersect(GrayVol_14(:,1),GrayVol_19(:,1));
  ID_inter=intersect(ThickAvg_14(:,1),ThickAvg_19(:,1));
 %ID_inter=intersect(SubcortexVol_14(:,1),SubcortexVol_19(:,1));


ID_inter=intersect(ID_inter,ID_506);

T1_14_inter=cd_align_intersect(ThickAvg_14,ID_inter);
T1_19_inter=cd_align_intersect(ThickAvg_19,ID_inter);


%% load TIV
ICV_14_inter=cd_align_intersect(ICV_14,ID_inter);
ICV_19_inter=cd_align_intersect(ICV_19,ID_inter);



%% cov_no_methylation_468
SDQ14=xlsread('./data/sdq14_ori_cov.xlsx');
cov_gender_site=[SDQ14(:,1),SDQ14(:,4:end)];

cov_no_methylation_468=cd_align_intersect(cov_gender_site,ID_inter);

cov_no_methylation_TIV_14=[cov_no_methylation_468,ICV_14_inter(:,2)];
cov_no_methylation_TIV_19=[cov_no_methylation_468,ICV_19_inter(:,2)];
cov_no_methylation_TIV=[cov_no_methylation_TIV_19;cov_no_methylation_TIV_14];


%% test thickavg
T1_19_14=[T1_19_inter;T1_14_inter];
grouplabel=cd_grouplabel(468,468);


[t,p,df]=cd_ttest2(T1_19_14(:,2:end), grouplabel, cov_no_methylation_TIV(:,2:end));
P_FDR=cd_FDR(p)




%%%%%%%%%%%%%%%%%%%%%%%
%% subcorex

 ID_inter=intersect(SubcortexVol_14(:,1),SubcortexVol_19(:,1));


ID_inter=intersect(ID_inter,ID_506);

T1_14_inter=cd_align_intersect(SubcortexVol_14,ID_inter);
T1_19_inter=cd_align_intersect(SubcortexVol_19,ID_inter);


%% load TIV
ICV_14_inter=cd_align_intersect(ICV_14,ID_inter);
ICV_19_inter=cd_align_intersect(ICV_19,ID_inter);



%% cov_no_methylation_468
SDQ14=xlsread('./data/sdq14_ori_cov.xlsx');
cov_gender_site=[SDQ14(:,1),SDQ14(:,4:end)];

cov_no_methylation_468=cd_align_intersect(cov_gender_site,ID_inter);

cov_no_methylation_TIV_14=[cov_no_methylation_468,ICV_14_inter(:,2)];
cov_no_methylation_TIV_19=[cov_no_methylation_468,ICV_19_inter(:,2)];
cov_no_methylation_TIV=[cov_no_methylation_TIV_19;cov_no_methylation_TIV_14];


T1_19_14=[T1_19_inter;T1_14_inter];
T1_19_14(:,1)=[];% dele ID


% subcortex order
ID_ROI_16=[15;13;6;12;8;7;5;1;29;28;24;27;26;25;23;19];% order by plot enigma

T1_19_14_ROI16=[];
for i=1:length(ID_ROI_16)
n=ID_ROI_16(i)

lin=T1_19_14(:,n);

T1_19_14_ROI16=[T1_19_14_ROI16,lin];

end



%% ttest
grouplabel=cd_grouplabel(468,468);


[t,p,df]=cd_ttest2(T1_19_14_ROI16, grouplabel, cov_no_methylation_TIV(:,2:end));
P_FDR=cd_FDR(p)

%% finished

close all; % close all figures
cd(working_dir);




%%











