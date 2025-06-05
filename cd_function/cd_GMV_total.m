function GMV_total=cd_GMV_total(dir)

% dir='/home2/migrate_handan/IMAGEN_20180305/IMAGEN_T1_Segment/T1_FU2_Raw_cat12/000000397377/mri/p1ADNIMPRAGEs003a1001.nii';

%% Method-1: Dpabi 

[d, VoxelSize, FileList, Header]=y_ReadAll(dir);
GMV_total=sum(d(:))*VoxelSize(1)*VoxelSize(2)*VoxelSize(3)/1000;

%% Method-2: load_untouch_nii ( from luojunyi)

% b=load_untouch_nii(dir);
% b.hdr.dime.dim
% b.hdr.dime.pixdim
% lin_size=b.hdr.dime.pixdim(2)*b.hdr.dime.pixdim(3)*b.hdr.dime.pixdim(4);
% img_mat=im2double(b.img);
% GMV_total=(sum(img_mat(:))*lin_size)/1000

end