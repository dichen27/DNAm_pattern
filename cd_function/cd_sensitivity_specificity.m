function  [accuracy,sensitivity,specificity,FPR]= cd_sensitivity_specificity(predict_label, ground_truth)  
%pat=1;control=0

accuracy=cd_accuracy(predict_label, ground_truth);

TP=0;
FN=0;
FP=0;
TN=0;
for i=1:length(predict_label)

%
if predict_label(i)==1 & ground_truth(i)==1
TP=TP+1;
else if predict_label(i)==0 & ground_truth(i)==1
FN=FN+1;
    else if predict_label(i)==1 & ground_truth(i)==0
FP=FP+1;
        else if predict_label(i)==0 & ground_truth(i)==0
TN=TN+1;
            end
        end
    end
end
end

%%
TPR=TP/(TP+FN); % Sensitivity
FPR=FP/(FP+TN); % 1-Specificity
TNR=TN/(TN+FP); % Specificity


sensitivity=TPR;
specificity=TNR;
end