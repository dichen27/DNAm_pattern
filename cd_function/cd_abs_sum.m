function norm1_someone=cd_abs_sum(bb)
% size of bb:  94 * 94

norm1_someone=[];
for j=1:size(bb,1)

norm_lin=bb(j,:);% some ROI

sum_lin=0;% for sum
for k=1:length(norm_lin)
    sum_lin=sum_lin+abs(norm_lin(k));   
end
mean_lin=sum_lin/(length(norm_lin)-1);% sum/n-1

norm1_someone=[norm1_someone,mean_lin];
end


end