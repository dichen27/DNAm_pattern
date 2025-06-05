%% svm_b
run_rate=[];
for j=0:3
for k=0:3
try
cmd=['-s ',num2str(j),' -t ',num2str(k)];
predict_final_b=[];
for i=1:size([fc_clu_b_zscore;fc_nor_zscore],1)
svm_fc_nor_b=p_t_nor_b_sort(1:4,3);% p<0.1

train_data_all_b=[fc_clu_b_zscore;fc_nor_zscore];
train_data_all_b(i,:)=[];% delete i people data
train_data_b=train_data_all_b(:,svm_fc_nor_b);

test_data_all_b=[fc_clu_b_zscore;fc_nor_zscore];
test_data_all_i=test_data_all_b(i,:);% select i people data
test_data_b=test_data_all_i(:,svm_fc_nor_b);

label_all=[]; %%%%算出trainlabel
for i=1:size(fc_clu_b_zscore,1)
    label_all=[label_all;-1];
end
for i=1:size(fc_nor_zscore,1)
    label_all=[label_all;1];
end%写出trainlabel

trainlabel_b=label_all;
trainlabel_b(i)=[];%删除第i个人的label
testlabel_b=label_all(i);%第i个人的label
model = svmtrain(trainlabel_b,train_data_b,cmd);
[predicted_label]=svmpredict(testlabel_b,test_data_b,model);%计算predicted_label
predict_final_b=[predict_final_b;predicted_label];
end


label_predict=[label_all,predict_final_b];
rate_num=0;
for i=1:length(label_all)
if label_predict(i,2)==label_predict(i,1)
    rate_num=rate_num+1;
end
end   
rate_b=rate_num/length(label_all);

run_rate_l={cmd,rate_b};
run_rate=[run_rate;run_rate_l]
catch
    continue;
end
end
end
% optimal select