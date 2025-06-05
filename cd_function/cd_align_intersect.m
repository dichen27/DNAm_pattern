function SDQ_14_intersect=cd_align_intersect(SDQ_14,ID_inter)


SDQ_14_intersect=[];
for i=1:length(ID_inter)
    index=SDQ_14(:,1)==ID_inter(i);
    lin=SDQ_14(index,:);
    
SDQ_14_intersect=[SDQ_14_intersect;lin];
end



end