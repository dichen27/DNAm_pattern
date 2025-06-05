

function m_zero=cd_cov_station(station_list)
Z_3=station_list;

max(Z_3);

% length(unique(score_pat12(:,9)))
m_zero=zeros(length(Z_3),max(Z_3));
% unique(head_clu_a(:,5))
for i=1:length(Z_3)
lie=Z_3(i);
m_zero(i,lie)=1;
end
m_zero(:,find(sum(abs(m_zero),1)==0))=[];%删除全是0的列
m_zero(:,size(m_zero,2))=[];%删除最后一列