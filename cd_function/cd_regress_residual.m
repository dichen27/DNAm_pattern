function residual_regressed=cd_regress_residual(Y,cov)
% Y=Drinking_14_intersect(:,2:end)
cov(:,find(sum(cov)==0))=[];% dele the all 0 cov_a_nor in cov

cov_one=[ones(size(cov,1),1),cov];

residual_regressed=[];
for i=1:size(Y,2)
Y_l=Y(:,i);
[b,bint,r,rint,stats] = regress(Y_l,cov_one);
residual_regressed=[residual_regressed,r];
end


end