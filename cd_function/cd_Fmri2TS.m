function  [TS,meanFD,ID]= cd_Fmri2TS(data_dir,mask_dir)  

% data_dir='C:/softc/abide2_test/';
% mask_dir='C:/softc/mat_tool/ABIDE/aal2_4mm.nii';
cd(data_dir)

sublist=dir();
ind=[sublist(:).isdir];
sublist=sublist(ind);
sublist=sublist(3:end);
n_sub=size(sublist,1);
disp(['Sample size of this site: ',num2str(n_sub),'.']);



TS={};
meanFD=[];
ID=[];
for k=1:length(sublist)
% k/length(sublist)

% Fmri to TS
cmd=sublist(k).name
fmri_cmd_dir=[data_dir,cmd,'/7_FunImg_to_Std/FunImg_4mmStdSpace_NoGlobalSignal.nii.gz'];
ts_lin=gong_fmri_ts(fmri_cmd_dir,mask_dir);
TS{k}=ts_lin;

% mean_FD
meanFD_dir=[data_dir,cmd,'/3_Motion_Corrected/meanFD_power.1D'];
meanFD_lin=textread(meanFD_dir);
meanFD=[meanFD;meanFD_lin];

% id
id_double_lin=str2num(cmd);
ID=[ID;id_double_lin];
end

end