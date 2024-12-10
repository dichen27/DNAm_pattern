function fc_new=cd_cell_ts_fc(data,r)
% nsub=length(ROI_ts);
% nroi=size(ROI_ts{1},2);
% nfc=(nroi-1)*nroi/2;
% 
% %network
% network=zeros(nsub,nfc);
% for i=1:nsub
%     net=corr(ROI_ts{i});
%     net=triu(net,1);
%     net(net==0)=[];
%     network(i,:)=net;
% end
% 
% network=0.5*log((1+network)./(1-network));

fc=[];
fc_new=[];

% 
% 
% data=ROI_12;
% r=94
bian=(r-1)*r/2;


fc=zeros(size(data,2),bian); 
for i=1:size(data,2)
aa=data{i};
aa=aa(:,1:r);
bb=corr(aa);
bb=bb-diag(diag(bb));%对角线变为0
kk=squareform(bb,'tovector');

fc(i,:)=kk;
end
fc_new=0.5.*log((1+fc)./(1-fc));%fisher变换