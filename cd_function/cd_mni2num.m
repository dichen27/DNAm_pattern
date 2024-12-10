function index_num=cd_mni2num(MNI_152,mask_dir)

% MNI_152=[-53,-33,44]
% mask_dir='/share/inspurStorage/home1/chendi/Documents/Depression/save_nii/Sub_hippo_t_map.nii';
% % i=2
% MNI_152=SampleCoordinates(i,2:4);
% mask_dir='/share/inspurStorage/home1/chendi/Documents/Depression/mask/Reslice_AAL3v1_121x145x121.nii';


mask_Img=load_nii(mask_dir);
T=[mask_Img.hdr.hist.srow_x;mask_Img.hdr.hist.srow_y;mask_Img.hdr.hist.srow_z;0,0,0,1];
index_num = cor_scan2vox_matlab(MNI_152,T);





% % try
% cor=MNI_152
% 
% % caution: if T is not given, the default T is
% if nargin == 1
%    %T = [-2,0,0,90;0,2,0,-126;0,0,2,-72;0,0,0,1];
%    error('Please provide the transformation matrix (T) .......')
% end
% 
% %cor = round(cor);
% Vox = inv(T)*[cor(:,1) cor(:,2) cor(:,3) ones(size(cor,1),1)]';
% Vox = Vox';
% Vox(:,4) = []
% % Vox=round(Vox)+1
% floor(Vox)+1

end

% cor=MNI_152

function Vox = cor_scan2vox_matlab(cor, T)
% function Vox = cor_scan2vox_matlab(cor, T)
% convert matrix coordinate to Vox index
%
% cor: an Nx3 matrix
% T: rotation matrix
% Vox is the returned coordinate in mni space
%
% caution: if T is not given, the default T is
if nargin == 1
   %T = [-2,0,0,90;0,2,0,-126;0,0,2,-72;0,0,0,1];
   error('Please provide the transformation matrix (T) .......')
end

%cor = round(cor);
Vox = inv(T)*[cor(:,1) cor(:,2) cor(:,3) ones(size(cor,1),1)]';
Vox = Vox';
Vox(:,4) = [];
Vox=round(Vox)+1;
% Vox=floor(Vox)+1;

end

