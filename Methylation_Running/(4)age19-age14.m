clear;
clc;
close;


addpath(genpath('/Users/dichen/Documents/MATLAB/cd_function/'));
working_dir='/Users/dichen/Documents/Methylation_Running/';
cd(working_dir);


%% load data
% data_dir='/share/inspurStorage/home1/chendi/Documents/Sylvane_Figure/methylation/';
% 
% t14=load('/home1/chendi/Documents/Sylvane_Figure/methylation/t14.txt');
% t_19=load([data_dir,'t_19.mat']);
% t19=t_19.t19;
% ID_506=t_19.t19_ID;
% 
% 
% % new working_dir methylation age
% working_dir='/share/inspurStorage/home1/chendi/Documents/Methylation_pattern/';
% cd(working_dir);

%% load data
data_dir='/Users/dichen/Documents/Methylation_Running/data_regressed/';


load([data_dir,'t14_full_regress.mat']);
load([data_dir,'t19_full_regress.mat']);
t14=t14_full_regress(:,2:end);
t19=t19_full_regress(:,2:end);

ID_506=t14_full_regress(:,1);


%% for cluster direct
substract_19_14_all_CpGs=[ID_506,t19-t14];
save("./List_ME/substract_19_14_all_CpGs.mat","substract_19_14_all_CpGs")

%r_substract_all_CpGs=corr(substract_19_14);
%save_dir=[working_dir,'List_ME/r_substract_all_CpGs.csv'];
%cd_save_csv(save_dir,r_substract_all_CpGs)

%%%%%%%%%%%%%%%%%%%%
%% load List consensus 19 & 14

list_consensus=load([working_dir,'List_Of_Methylation_Change_Power5_Min50_Merge0_UnsignedTOM_Consensus_19_14.txt']);


substract_19_14=t19-t14;
num=unique(list_consensus);


ME_consensus_14=[];
ME_consensus_19=[];
ME_consensus_substract_19_14=[];
num_ME=[];
for i=1:length(num)
   
cmd=num(i)
    
index=list_consensus==cmd;
   

% age 14
data_14=t14(:,index);
lin_14=mean(data_14,2);
ME_consensus_14=[ME_consensus_14,lin_14];


% age 19
data_19=t19(:,index);
lin_19=mean(data_19,2);
ME_consensus_19=[ME_consensus_19,lin_19];

% age19-14
data_19_14=substract_19_14(:,index);
lin_19_14=mean(data_19_14,2);
ME_consensus_substract_19_14=[ME_consensus_substract_19_14,lin_19_14];



num_ME=[num_ME;sum(index)];
end
ME_consensus_14(:,1)=[];% dele ME-0
ME_consensus_19(:,1)=[];
ME_consensus_substract_19_14(:,1)=[];

%%


r_14=corr(ME_consensus_14);
r_19=corr(ME_consensus_19);

r_14_19=corr(ME_consensus_14,ME_consensus_19);

r_substract=corr(ME_consensus_substract_19_14);


save_dir=[working_dir,'List_ME/ME_consensus_14.csv'];
cd_save_csv(save_dir,ME_consensus_14)


save_dir=[working_dir,'List_ME/ME_consensus_19.csv'];
cd_save_csv(save_dir,ME_consensus_19)

save_dir=[working_dir,'List_ME/ME_consensus_substract_19_14.csv'];
cd_save_csv(save_dir,ME_consensus_substract_19_14)

save_dir=[working_dir,'List_ME/r_14.csv'];
cd_save_csv(save_dir,r_14)


save_dir=[working_dir,'List_ME/r_19.csv'];
cd_save_csv(save_dir,r_19)

save_dir=[working_dir,'List_ME/r_14_19.csv'];
cd_save_csv(save_dir,r_14_19)

save_dir=[working_dir,'List_ME/r_substract.csv'];
cd_save_csv(save_dir,r_substract)





%% load index_450k_850k
load([working_dir,'index_450k_850k.mat']);
sum(index_450k)
sum(index_850k)

% age 14
ME_consensus_14_450k=ME_consensus_14(index_450k,:);
ME_consensus_14_850k=ME_consensus_14(index_850k,:);


save_dir=[working_dir,'List_ME/ME_consensus_14_450k.csv'];
cd_save_csv(save_dir,ME_consensus_14_450k)
save_dir=[working_dir,'List_ME/ME_consensus_14_850k.csv'];
cd_save_csv(save_dir,ME_consensus_14_850k)

% age 19
ME_consensus_19_450k=ME_consensus_19(index_450k,:);
ME_consensus_19_850k=ME_consensus_19(index_850k,:);


save_dir=[working_dir,'List_ME/ME_consensus_19_450k.csv'];
cd_save_csv(save_dir,ME_consensus_19_450k)

save_dir=[working_dir,'List_ME/ME_consensus_19_850k.csv'];
cd_save_csv(save_dir,ME_consensus_19_850k)

% age19 - age14
ME_consensus_substract_19_14_450k=ME_consensus_substract_19_14(index_450k,:);
ME_consensus_substract_19_14_850k=ME_consensus_substract_19_14(index_850k,:);

save_dir=[working_dir,'List_ME/ME_consensus_substract_19_14_450k.csv'];
cd_save_csv(save_dir,ME_consensus_substract_19_14_450k)
save_dir=[working_dir,'List_ME/ME_consensus_substract_19_14_850k.csv'];
cd_save_csv(save_dir,ME_consensus_substract_19_14_850k)


%% finished
cd(working_dir);


















