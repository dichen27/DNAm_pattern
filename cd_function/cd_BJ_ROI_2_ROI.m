function [ID_baseline_BJ,FC_seed]=cd_BJ_ROI_2_ROI(data_baseline_dir,mm_ROI_baseline,mm_ROI_FU)

%% ROI to ROI

% data_baseline_dir='/home1/chendi/Cao_ADHD/data/';
cd(data_baseline_dir);
sublist=dir(data_baseline_dir);
sublist(1:2,:)=[];


FC_seed=[];
ID_baseline_BJ=[];
for k=1:length(sublist)

k
cmd=sublist(k).name;
% disp(['Name of ID: ',cmd]);

data_str=load_nii([data_baseline_dir,cmd,'/Filtered_4DVolume.nii']);


% data_str=load_nii('F:\ABCD_1\NDARINVBDMGZK3B_baselineYear1Arm1_ABCD-MPROC-rsfMRI_20170615165443\7_FunImg_to_Std\FunImg_3mmStdSpace_NoGlobalSignal.nii.gz');
img_lin=data_str.img;


% TS ROI_baseline
ts_ROI_baseline=[];
for i=1:size(img_lin,4)

img_ti=data_str.img(:,:,:,i);
img_map_ti=double(img_ti(mm_ROI_baseline)');

mean_lin=nanmean(img_map_ti);
ts_ROI_baseline=[ts_ROI_baseline;mean_lin];
end


% TS ROI_FU
ts_ROI_FU=[];
for i=1:size(img_lin,4)

img_ti=data_str.img(:,:,:,i);
img_map_ti=double(img_ti(mm_ROI_FU)');

mean_lin=nanmean(img_map_ti);
ts_ROI_FU=[ts_ROI_FU;mean_lin];
end


% TS to FC

kk=corr(ts_ROI_baseline,ts_ROI_FU);

kk_fisher=0.5.*log((1+kk)./(1-kk));%fisher�z-form



FC_seed=[FC_seed;kk_fisher];
ID_baseline_BJ=[ID_baseline_BJ;cmd];
end







end