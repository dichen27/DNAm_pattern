




%% LDA
% input_sample=(testdata_all(:,a_u)-repmat(mean(testdata_all(:,a_u)),1,1))*c;%计算testdata
% input_sample=testdata_all(:,a_u);
input_sample=U_test(:,1);

% training=(train_data_all(:,a_u)-repmat(mean(train_data_all(:,a_u)),size(train_data_all,1),1))*c;
% training=train_data_all(:,a_u);
training=U_train(:,1);


train_clu=[training,clu];

index=train_clu(:,size(train_clu,2))==1;
class1_samples=train_clu(index,1:(size(train_clu,2)-1));
index=train_clu(:,size(train_clu,2))==2;
class2_samples=train_clu(index,1:(size(train_clu,2)-1));

[H,p,ci,stats]=ttest2(class1_samples,class2_samples)

output_class=cd_LDA(input_sample, class1_samples, class2_samples);
tabu=tabulate(output_class)

THR=min(tabu(:,2))
label_dream=[(ones(THR,1));2*ones((length(filename_pat_cv)-THR),1)];
label_predict=[label_dream,w_cv(:,10)];
rate_num=0;
for i=1:length(label_dream)
if label_predict(i,2)==label_predict(i,1)
    rate_num=rate_num+1;
end
end   
rate=rate_num/length(label_dream)
%% SVM
pval=[];
rho=[];
rate=[];
for dong=0:3
%    dong=0
cmd=['-s 0 -t ',num2str(dong)];
% cmd=['-s 0 -t ',num2str(t)];

% train_data=train_data_all(:,a_u);
train_data=U_train;

train_label=clu;


% test_data=testdata_all(:,a_u);
test_data=U_test;

test_label=[];
for i=1:length(filename_pat_cv)
if score_double12(i,6)>14
    test_label_l=1;
else
    test_label_l=2;
end
test_label=[test_label;test_label_l];
end


model = svmtrain(train_label,train_data,cmd);
[predicted_label]=svmpredict(test_label,test_data,model);%计算predicted_label

score_double12_label=[score_double12,predicted_label];
[rho_l,pval_l]=corr(score_double12_label(:,6),score_double12_label(:,10));
pval=[pval;pval_l];
rho=[rho;rho_l];
w_cv=sortrows(score_double12_label,-6);
end
rate
pval
rho

%% Classify
training=(train_data_all(:,a_u)-repmat(mean(train_data_all(:,a_u)),size(train_data_all,1),1))*c;

group=clu;
testdata_all=fc_pat_test;
% U=testdata_all(:,a_u);
sample=(testdata_all(:,a_u)-repmat(mean(testdata_all(:,a_u)),1,1))*c;%计算testdata

class = classify(sample,training,group);
tabu=tabulate(class)


score_double12_label=[score_double12,class];
[rho_l,pval_l]=corr(score_double12_label(:,6),score_double12_label(:,10))
w_cv=sortrows(score_double12_label,-6);


THR=min(tabu(:,2))
label_dream=[(ones(THR,1));2*ones((length(filename_pat_cv)-THR),1)];
label_predict=[label_dream,w_cv(:,10)];
rate_num=0;
for i=1:length(label_dream)
if label_predict(i,2)==label_predict(i,1)
    rate_num=rate_num+1;
end
end   
rate=rate_num/length(label_dream)
%% linkage
testdata_all=fc_pat_test;
% U=testdata_all(:,a_u);
U=(testdata_all(:,a_u)-repmat(mean(testdata_all(:,a_u)),1,1))*c;%计算testdata

U_v=pdist(U);
% figure;histogram(U_v);
U_m=squareform(U_v);
U_lin=linkage(U_m,'ward');
figure;dendrogram(U_lin,0);%层次聚类

testn=[];%选择聚类个数
for i=1:10    
clu_t=cluster(U_lin,i);
s=silhouette(U,clu_t);
testn=[testn;mean(s)];
end
testn
testn_max=max(testn);
[testn_max_e,testn_f]=find(testn==testn_max);
clu_linkage=cluster(U_lin,2);%具体看聚类结果
tabu=tabulate(clu_linkage)
%  kankan=[ scatter_c,clu];

score_double12_label=[score_double12,clu_linkage];
[rho_l,pval_l]=corr(score_double12_label(:,6),score_double12_label(:,10))
w_cv=sortrows(score_double12_label,-6);

THR=min(tabu(:,2))
label_dream=[(ones(THR,1));2*ones((length(filename_pat_cv)-THR),1)];
label_predict=[label_dream,w_cv(:,10)];
rate_num=0;
for i=1:length(label_dream)
if label_predict(i,2)==label_predict(i,1)
    rate_num=rate_num+1;
end
end   
rate=rate_num/length(label_dream)
%% clusterdata
testdata_all=fc_pat_test;
% test_data=testdata_all(:,a_u);
test_data=(testdata_all(:,a_u)-repmat(mean(testdata_all(:,a_u)),1,1))*c;%计算testdata

clu_cluster= clusterdata(test_data,'linkage','ward','savememory','on','maxclust',2);
 figure;scatter3(test_data(:,1),test_data(:,2),test_data(:,3),20,clu_cluster,'filled')
 

 
tabu=tabulate(clu_cluster)

score_double12_label=[score_double12,clu_cluster];
[rho_l,pval_l]=corr(score_double12_label(:,6),score_double12_label(:,10))
w_cv=sortrows(score_double12_label,-6);


THR=min(tabu(:,2))
% THR=23
label_dream=[(ones(THR,1));2*ones((length(filename_pat_cv)-THR),1)];
label_predict=[label_dream,w_cv(:,10)];
rate_num=0;
for i=1:length(label_dream)
if label_predict(i,2)==label_predict(i,1)
    rate_num=rate_num+1;
end
end   
rate=rate_num/length(label_dream)

