function [PPR,PPV,NPV]=cd_PPR(predict_label,groud_truth)



% predict_label=cd_grouplabel(length(DAWBA_pat14),length(DAWBA_control_14));
% groud_truth=[DAWBA_pat14(:,2);DAWBA_control_14(:,2)];


data=[predict_label,groud_truth];
data_sort=sortrows(data,-1);


index=data_sort(:,1)==1;
pat=data_sort(index,2);

index=data_sort(:,1)==0;
control=data_sort(index,2);

PPV=sum(pat)/length(pat);
NPV=1-sum(control)/length(control);
PPR=PPV/(1-NPV);


% PPV_14=sum(DAWBA_pat14(:,2))/size(DAWBA_pat14,1)
% NPV_14=1-sum(DAWBA_control_14(:,2))/size(DAWBA_control_14,1)
% 
% PPR_14=PPV_14/(1-NPV_14)

end