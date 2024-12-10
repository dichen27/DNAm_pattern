function img=cd_img2dto3dmask(size3dimage,dim_vector)

load('/share/home1/chendi/Documents/14MID_large_no/data14_activate');

img=zeros(53,63,46);
for i=1:size(data14_activate,1)
    index=sample==data14_activate(i,4);
img(index)=1;
%     img(dim_vector(i,1),dim_vector(i,2),dim_vector(i,3))=1;
end
return