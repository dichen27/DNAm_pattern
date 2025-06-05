function ROI_ts=g_extract_ROI_ts(fMRI_path,mask_path)
%        ROI_ts=g_extract_ROI_ts(fMRI_path,mask_path)


a=load_nii(mask_path);
mask1=a.img;
roi1=unique(mask1);
roi1(roi1==0)=[];
nroi=length(roi1);

aa=load_nii(fMRI_path);
image=aa.img;
nts=size(image,4);


ROI_ts=zeros(nts,nroi);
mm=(mask1>=min(roi1) & mask1<=max(roi1));
mm1=mask1(mm);
image2d=zeros(sum(mm(:)),nts);
for j=1:nts
    img=image(:,:,:,j);
    image2d(:,j)=img(mm);
end
for j=1:nroi
    ROI_ts(:,j)=mean(image2d(mm1==roi1(j),:),'omitnan');
end


end