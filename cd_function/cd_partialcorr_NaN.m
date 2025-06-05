
function [r,p,df,n]=cd_partialcorr_NaN(PRS_mean_ID,SDQ_total_4groups,cov_GWAS)

ID_intersect=intersect(PRS_mean_ID(:,1),SDQ_total_4groups(:,1));
ID_intersect=intersect(ID_intersect,cov_GWAS(:,1));

PRS_mean_intersect=[];
for i=1:length(ID_intersect)

index=PRS_mean_ID(:,1)==ID_intersect(i);
    lin=PRS_mean_ID(index,:);
    
  PRS_mean_intersect=[PRS_mean_intersect;lin];  
end


SDQ_intersect=[];
for i=1:length(ID_intersect)
index=SDQ_total_4groups(:,1)==ID_intersect(i);
    lin=SDQ_total_4groups(index,:);

    SDQ_intersect=[SDQ_intersect;lin];
end

cov_intersect=[];
for i=1:length(ID_intersect)
index=cov_GWAS(:,1)==ID_intersect(i);
    lin=cov_GWAS(index,:);

    cov_intersect=[cov_intersect;lin];
end


df=length(ID_intersect)-2-size(cov_intersect,2)+1;

n=length(ID_intersect);
[r,p]=partialcorr(PRS_mean_intersect(:,2:end),SDQ_intersect(:,2:end),cov_intersect(:,2:end),'rows','complete');

end