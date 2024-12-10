function [log_p_all]=cd_log_p(p_all)



log_p_all=[];
for j=1:size(p_all,2)

log_p_all_l=[];
for i=1:length(p_all(:,j))
log_p_all_ll=(-1)*log10(p_all(i,j));
log_p_all_l=[log_p_all_l;log_p_all_ll];
end

log_p_all=[log_p_all,log_p_all_l];
end





end