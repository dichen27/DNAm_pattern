function [block]=cd_find_block(j,FC)%iΪ������С��3Ϊ120С����4Ϊ��д,3Ϊ������

% [~,~,AAL2] =xlsread('AAL2.xlsx');

load('aal2_name.mat');
AAL2=aal2_name;

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

symbol=char();%symbol'--'
symbol='��';

block={};
% p=char();%���a������
% q=char()
for i=1:size(FC);
number_a=a(i);
p_t=AAL2{number_a,j};
%a_t_br=AAL2{number,3};
%a_tt=[a_t,a_t_br]
number_b=b(i);
q_t=AAL2{number_b,j};
block_t=[p_t,symbol,q_t];
block=[block;block_t];
end

end
