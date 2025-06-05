
function [sample_120]=cd_map_tval(t_location,num)
[~,~,AAL2] =xlsread('AAL2.xlsx');



id_matrix=squareform(1:4371);
sample=zeros(94,94);


for i=1:num

sample(find(id_matrix==t_location(i,2)))=t_location(i,1);

end
sample_120=zeros(120,120);
sample_120(1:94,1:94)=sample;
end