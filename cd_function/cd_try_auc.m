function [result]=cd_try_auc(output,test_targets)

%计算AUC值,test_targets为原始样本标签,output为分类器得到的标签
%均为行或列向量

% test_targets=test_label_clinic_clu;
% output=predict_label;

[A,I]=sort(output);
M=0;N=0;
for i=1:length(output)
if(test_targets(i)==1)
M=M+1;
else
N=N+1;
 end
end
sigma=0;
for i=M+N:-1:1
 if(test_targets(I(i))==1)
 sigma=sigma+i;
end
end
result=(sigma-(M+1)*M/2)/(M*N);



% 
% --------------------- 
% 作者：Roc-Ng 
% 来源：CSDN 
% 原文：https://blog.csdn.net/windows_peng/article/details/70833158 
% 版权声明：本文为博主原创文章，转载请附上博文链接！