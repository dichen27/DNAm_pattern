function [result]=cd_try_auc(output,test_targets)

%����AUCֵ,test_targetsΪԭʼ������ǩ,outputΪ�������õ��ı�ǩ
%��Ϊ�л�������

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
% ���ߣ�Roc-Ng 
% ��Դ��CSDN 
% ԭ�ģ�https://blog.csdn.net/windows_peng/article/details/70833158 
% ��Ȩ����������Ϊ����ԭ�����£�ת���븽�ϲ������ӣ�