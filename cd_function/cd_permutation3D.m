function [positive_cluster,negative_cluster,ts,mask1,mask2]=permutation3D(X,design,mask,CDT,nperm)
%        [positive_cluster,negative_cluster,ts,mask1,mask2]=permutation3D(X,design,mask,CDT,nperm)
%Input:
%      X: n-by-m matrix, n sample and m voxels
%      design: n-by-p matrix, n sample and p covariates (The first column
%              is the target variable), do not need to include a column of
%              1 in the design matrix.
%      mask: non-zero elements will be used
%      CDT: cluster defining threshold p-value (default=0.01), allow multiple
%           thresholds in a vector(e.g. [0.05,0.01,0.005,0.001])
%      nperm: number of permutation (default=1000)
%Output:
%      positive_cluster: the cluster size of positively correlated clusters and their
%                         p-values
%      negative_cluster: the cluster size of negatively correlated clusters and their
%                         p-values
%
% Written by Weikang Gong (E-mail: weikanggong@gmail.com)
%
if nargin<4
    CDT=0.01;
    nperm=1000;
end
 
if nargin<5
    nperm=1000;
end
 
[a1,a2,a3]=size(mask);
[nn,~]=size(X);
qq=size(design,2);%????????????
[d1,d2,d3]=ind2sub(size(mask),find(mask~=0));
dim=[d1,d2,d3];
mm=length(CDT);%????CDT??????
thre1=abs(tinv(CDT,nn-qq-1));
thre2=tinv(CDT,nn-qq-1);
%Permutation SPM
maxcl1=zeros(nperm,mm);
maxcl2=zeros(nperm,mm);
for i=1:nperm
   
    ind=datasample(1:nn,nn,'Replace',false);
    ts=BWAS_Tregression(design(ind,:),X);
    for j=1:mm
        ind1=ts>thre1(j);
        ind2=ts<thre2(j);
        mask1=img2dto3dmask([a1,a2,a3],dim(ind1,:));
        mask2=img2dto3dmask([a1,a2,a3],dim(ind2,:));
        cc1=bwconncomp(mask1,18);
        cc2=bwconncomp(mask2,18);
        ss1=max(cellfun(@length,cc1.PixelIdxList));
        ss2=max(cellfun(@length,cc2.PixelIdxList));
        if isempty(ss1)
            maxcl1(i,j)=0;
        else
            maxcl1(i,j)=ss1;
        end
        if isempty(ss2)
            maxcl2(i,j)=0;
        else
            maxcl2(i,j)=ss2;
        end
    end
    i
end
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%??????????????????permutatiopn????????maxcl1??maxcl2????????????????????????????????cluster????????????????????
%??????????????cluster????????????????????????????cluster??????????????p-value
ts=BWAS_Tregression(design,X);   %ts????????t??????
positive_cluster=cell([mm,1]); 
negative_cluster=cell([mm,1]);
for j=1:mm
    positive_cluster{j,2}=thre1(j);
    negative_cluster{j,2}=thre2(j);
    positive_cluster{j,3}=CDT(j);
    negative_cluster{j,3}=CDT(j);
    ind1=ts>thre1(j);
    ind2=ts<thre2(j);
    mask1=img2dto3dmask([a1,a2,a3],dim(ind1,:));
    mask2=img2dto3dmask([a1,a2,a3],dim(ind2,:));
    cc1=bwconncomp(mask1,18);
    cc2=bwconncomp(mask2,18);
    ss1=sort(cellfun(@length,cc1.PixelIdxList),'descend');
    ss2=sort(cellfun(@length,cc2.PixelIdxList),'descend');
    if ~isempty(ss1)
        for jj=1:length(ss1)
            positive_cluster{j,1}(jj,1)=ss1(jj);
            positive_cluster{j,1}(jj,2)=mean(maxcl1(:,j)>ss1(jj));%????????????p-value
        end
    else
        positive_cluster{j,1}=[];
    end
    if ~isempty(ss2)
        for jj=1:length(ss2)
            negative_cluster{j,1}(jj,1)=ss2(jj);
            negative_cluster{j,1}(jj,2)=mean(maxcl2(:,j)>ss2(jj));
        end
    else
        negative_cluster{j,1}=[];
    end
end
  
return
 

 
 
 