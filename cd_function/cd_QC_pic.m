
function  [TS,meanFD,ID]= cd_QC_pic(data_dir,output_dir)  


% data_dir=data_dir_lin
% output_dir='/home1/chendi/Documents/qc_newdata_pic/';
cd(data_dir);


sublist=dir();
ind=[sublist(:).isdir];
sublist=sublist(ind);
sublist=sublist(3:end);
n_sub=size(sublist,1);
disp(['Sample size of this site: ',num2str(n_sub),'.']);



% copy to output_dir

for k=1:length(sublist)

cmd=sublist(k).name

pic_dir_someone=[data_dir,cmd,'/',cmd,'_normalization.tif'];
unix(['cp ',pic_dir_someone,' ',output_dir]);

end


end

% pic_dir_someone='/home1/chendi/Documents/abide1/fMRI/Newdata_abide1/Caltech/51456/51456_normalization.tif';
% output_dir='/home1/chendi/Documents/qc_newdata_pic/';
% 
% 
% unix(['cp ',pic_dir_someone,' ',output_dir]);
