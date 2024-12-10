function [Covariance,intercept,W]=cd_WBC_covariance_improve(VDI,t_12,N_1,N_2,h_1,h_2,intercept_1,intercept_2,N_s)

% VDI=VID_MID_SST_T1;
% 
% t_12=t_12_all;
% N_s=min(N_1,N_2)

%%

M=length(VDI);
%% first regression

y=t_12;
x=VDI;



x_one=[ones(size(x,1),1),x];
[beta,bint,r,rint,stats] = regress(y,x_one);

beta_foundation=beta(2);

Q_g_foundation=(beta(2)*M)/sqrt(N_1*N_2)
intercept_foundation=beta(1)



%% weight for second regression

% coefficient
coefficient_1=(N_1*h_1)/M;
coefficient_2=(N_2*h_2)/M;
Q_g_lin=sum(t_12)/(N_1*N_2);

coefficient_3=(Q_g_lin*sqrt(N_1*N_2))/M;

W=zeros(M,M);
for i=1:length(VDI)
    
var_1=(coefficient_1*VDI(i)+intercept_1)*(coefficient_2*VDI(i)+intercept_2);

% if N_s==0;
%     var_2=(Q_g_foundation*sqrt(N_1*N_2)*VDI(i))/M;
% else
var_2=((Q_g_foundation*sqrt(N_1*N_2)*VDI(i))/M)+intercept_foundation;
var_2=beta(2)*VDI(i)+beta(1);
% end

var_all=var_1+(var_2)^2;

lin=(1/var_all)*(1/VDI(i));    %(1/var)*(1/L_j)
% lin=1/(R_2_moment(i)+1);
% lin=(1/var_all);
W(i,i)=lin;
end



%% second regression

x=VDI;
X=[ones(length(x),1),x];


lin=inv(X'*W*X);
beta=lin*X'*W*y;
% Res=y-X*beta_final;

beta_final=beta(2);

intercept=beta(1);

Covariance=(beta_final*M)/(sqrt(N_1*N_2));

%% Genetic_Correlation
% r=Covariance/sqrt(h_1*h_2);



end




