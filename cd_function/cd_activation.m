function data_14_mean_noregerss=cd_activation(brain)
% size of bb:  94 * 94
data_somegroup=brain;

data_14_mean_noregerss=[];% tval
for i=1:size(data_somegroup,2)

%    i/size(data_somegroup,2)
   
     
data_somevoxel=data_somegroup(:,i);

[h,p,ci,stats]=ttest(data_somevoxel);
ef_lin=stats.tstat/sqrt(stats.df);


%% one-sample ttest by myself
% [m_l,n_l]=find(isnan(data_somevoxel)==1);% rm NaN
% data_somevoxel(m_l,:)=[];
% 
% 
% 
% somevoxel_mean=mean(data_somevoxel);
% somevoxel_std=std(data_somevoxel);
% somevoxel_sqrt=sqrt(length(data_somevoxel)-1);% 
% 
% % df = max(length(data_somevoxel) - ones('like',data_somevoxel), 0); % the df from ttest function
% 
% somevoxel_14_mean=(somevoxel_mean*somevoxel_sqrt)/somevoxel_std;%(mean*sqrt)/std
%%

data_14_mean_noregerss=[data_14_mean_noregerss;ef_lin];

end




end