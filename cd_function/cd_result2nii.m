function cd_result2nii(ROI_result,mask_dir,save_dir)


% mask_dir='/share/inspurStorage/home1/chendi/Documents/Depression/mask/AAL3v1_1mm.nii';
% save_dir='/share/inspurStorage/home1/chendi/Documents/Depression/save_nii/ACC_stage/ACC_stage_1.nii';

% AAL mask
mask=load_nii(mask_dir);
mask_img=mask.img;
mask_img_double=double(mask_img);
roi1=unique(mask_img_double);
roi1(roi1==0)=[];
nroi=length(roi1);


mm=(mask_img_double>=min(roi1) &mask_img_double<=max(roi1));
block_label=mask_img_double(mm);% prepare block_label
ROI_tabulate=tabulate(block_label);



img_result=zeros(size(mask_img_double));
for i=1:size(ROI_result)

index=mask_img_double==ROI_result(i,1);
img_result(index)=ROI_result(i,2);
end


% save
% someone_mask=load_nii([working_dir,'mask/Reslice_AAL3v1_121x145x121.nii']);

mask.hdr.dime.datatype = 64;
mask.hdr.dime.bitpix = 64;
mask.filetype=64;

mask.img=img_result;



save_nii(mask,save_dir);

end







