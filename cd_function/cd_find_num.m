function [kkk]=cd_find_block(FC)%iΪ������С��1Ϊ120С����3Ϊ������
[~,~,AAL2] =xlsread('AAL2.xlsx');




id_matrix=squareform(1:4371);
id_matrix(find(triu(id_matrix,1)==0))=0;



a=[];
b=[];
for i=1:size(FC)
[e,f]=find(id_matrix==FC(i));
a=[a;e];
b=[b;f];
end
kkk=[a,b];
end
