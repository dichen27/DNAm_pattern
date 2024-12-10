function [t_DSM,Sample_size,score_index,p_all,r_all]=cd_VDI_t(brain_id,score_id,cov_id)




% Sample_size=size(score_id,1);

Sample_size=[];
score_index=[];
t_DSM=[];
p_all=[];
r_all=[];
for k=2:size(score_id,2)
    
%  k
t_k=[];
p_k=[];
r_k=[];
Sample_size_k=[];
for  i=2:size(brain_id,2)
   
data_lin=[brain_id(:,i),score_id(:,k),cov_id(:,2:end)];
data_lin(any(isnan(data_lin), 2),:) = [];% dele any row NaN

sample_size_lin=size(data_lin,1);
    
    
[r_lin,p_lin]=partialcorr(data_lin(:,1),data_lin(:,2),data_lin(:,3:end));

% r to t

n_lin=sample_size_lin-size(data_lin(:,3:end),2);% sample size-- cov
t_lin=cd_r2t(r_lin,n_lin);

Sample_size_k=[Sample_size_k;sample_size_lin];
t_k=[t_k;t_lin];
p_k=[p_k;p_lin];
r_k=[r_k;r_lin];
end
Sample_size_mean_k=mean(Sample_size_k);
t_DSM=[t_DSM,t_k];
score_index=[score_index;k];
Sample_size=[Sample_size;Sample_size_mean_k];
p_all=[p_all,p_k];
r_all=[r_all,r_k];
end

















end