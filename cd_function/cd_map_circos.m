
function [sample]=cd_map_circos(val_location)%iΪ������С��1Ϊ120С����3Ϊ������
[~,~,AAL2] =xlsread('AAL2.xlsx');



id_matrix=squareform(1:4371);
sample=zeros(94,94);


for i=1:size(val_location,1)

sample(find(id_matrix==val_location(i,2)))=val_location(i,1);

end
end