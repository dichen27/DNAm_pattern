clear;
clc;
close;

addpath(genpath('/Users/dichen/Documents/MATLAB/cd_function'));
working_dir='/Users/dichen/Documents/Methylation_Running/';
cd(working_dir);

%% load only for ID_506: ID_506 at age14 the same with age19
t_19=load('./data/t_19.mat');
t19=t_19.t19;
ID_506=t_19.t19_ID;


%% cov methylation
cov_methy_14=csvread('./data/methy_cov14_506.csv',1,0);
% cov_methy_14_inter(:,10:end)=[];



%% load DNAm-14
ME_14_all=csvread('./data_14_origin_cluster18.csv',1,1);
Module0_14=csvread('(Reply1-2)data_14_module_0.csv',1,1);
Module0_14_ID=[ID_506,Module0_14];

%% 10 brain-related clusters
select_cluster=[1,3:9,11,12]';
length(select_cluster)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% subtance use
output=[];
P_all=[];
r_all=[];
for n=1:length(select_cluster)
k=select_cluster(n);

ME_14_all_ID=[ID_506,ME_14_all(:,k)];


data_all=importdata('./BL_Behaviour/IMAGEN-IMGN_ESPAD_CHILD_RC5-IMAGEN_DIGEST.csv');
phenotype_14=data_all.data;
name_phenotype=data_all.textdata(2:end);

[r,p,df,n]=cd_partialcorr_NaN(ME_14_all_ID,phenotype_14,cov_methy_14);
lin=[r',p'];

output=[output,lin];
P_all=[P_all,p'];
r_all=[r_all,r'];
end

% FDR-correction
P_fdr=cd_FDR(P_all);


%% subtance use module-0
data_all=importdata('./BL_Behaviour/IMAGEN-IMGN_ESPAD_CHILD_RC5-IMAGEN_DIGEST.csv');
phenotype_14=data_all.data;
name_phenotype=data_all.textdata(2:end)

[r,p_substance,df,n]=cd_partialcorr_NaN(Module0_14_ID,phenotype_14,cov_methy_14)


%% MDD module-0
data_all=importdata('./IMAGEN_dawba_cleaned_BL.xlsx');
phenotype_14=data_all.data(:,1:2);

[r,p_MDD,df,n]=cd_partialcorr_NaN(Module0_14_ID,phenotype_14,cov_methy_14)
P_module0_14=[p_substance';p_MDD];
save('./(A-2)P_module0_14.mat','P_module0_14')



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MDD with DNAm-14
ME_14_all_ID=[ID_506,ME_14_all(:,select_cluster)];

data_all=importdata('./IMAGEN_dawba_cleaned_BL.xlsx');
phenotype_14=data_all.data(:,1:2);

[r,p,df,n]=cd_partialcorr_NaN(ME_14_all_ID,phenotype_14,cov_methy_14)
P_dfr=cd_FDR(p)

%% save for steiger
% intersect
ID_inter=intersect(ME_14_all_ID(:,1),phenotype_14(:,1));
ID_inter=intersect(ID_inter,cov_methy_14(:,1));

phenotype_14_inter=cd_align_intersect(phenotype_14,ID_inter);
ME_14_select_ID=cd_align_intersect(ME_14_all_ID,ID_inter);
cov_methy_14_inter=cd_align_intersect(cov_methy_14,ID_inter);


% residual
MDD_14_residulized=cd_regress_residual(phenotype_14_inter(:,2),cov_methy_14_inter(:,2:9));% sex,site
ME_14_cluster_7=cd_regress_residual(ME_14_select_ID(:,7),cov_methy_14_inter(:,2:end));



ME_14_residulized_cluster_7_ID=[ID_inter,ME_14_cluster_7];
MDD_14_residulized=[ID_inter,MDD_14_residulized];

%% save for steiger
save('./(A-2)save_for_steiger_MDD_cluster7_14.mat','ME_14_residulized_cluster_7_ID','MDD_14_residulized')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% finished
close all; % close all figures
cd(working_dir);














