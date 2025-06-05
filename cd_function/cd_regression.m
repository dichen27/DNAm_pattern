function [T,P,Res]=cd_regression(X,Y)


 
contrast=[1,zeros(1,size(X,2))]';
 
X=[X,ones(size(X,1),1)];
 
df=size(X,1)-size(X,2);
 
ss=inv(X'*X);
 
beta=ss*X'*Y;
 
Res=Y-X*beta;
 
sigma=sqrt(sum(Res.^2)./df);
 
T=beta'*contrast./(sigma'*sqrt(contrast'*ss*contrast));
 
%stdRes=Res./sigma;
 
P=2*(1-tcdf(abs(T),df));
 
 
return