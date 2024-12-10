function [Covariance,intercept]=cd_WBC_covariance(VDI,t_12,N_1,N_2,h_1,h_2,intercept_1,intercept_2,N_s)

%% prepare data




%% update algorithm
M=length(VDI);
%% first regression improve

% h_1=max(h_1,0);
% h_1=min(h_1,1);
% 
% 
% h_2=max(h_2,0);
% h_2=min(h_2,1);

% var part 1
coefficient_1=(N_1*h_1)/M;
coefficient_2=(N_2*h_2)/M;


Q_g_lin=sum(t_12)/(N_1*N_2);
Q_g_lin=max(Q_g_lin,-1);
Q_g_lin=min(Q_g_lin,1);

coefficient_3=(Q_g_lin*sqrt(N_1*N_2))/M;


% weight for first regression

W=zeros(M,M);
for i=1:length(VDI)
    
var_1=(coefficient_1*VDI(i)+intercept_1)*(coefficient_2*VDI(i)+intercept_2);
var_2=(coefficient_3*VDI(i))^2;
var_all=var_1+var_2;

% lin=1/var_all;    % 1/var

lin=(1/var_all)*(1/VDI(i));    % 1/L_j

W(i,i)=lin;
end



y=t_12;
x=VDI;
X=[ones(length(x),1),x];

lin=inv(X'*W*X);
beta=lin*X'*W*y;
% Res=y-X*beta;

beta_found=beta;

Q_g_foundation=(beta(2)*M)/sqrt(N_1*N_2);
% Q_g_foundation=max(Q_g_foundation,-1);
% Q_g_foundation=min(Q_g_foundation,1);

% weight for second regression
W=zeros(M,M);
for i=1:length(VDI)
var_1=(coefficient_1*VDI(i)+intercept_1)*(coefficient_2*VDI(i)+intercept_2);

if N_s==0;
    var_2=(Q_g_foundation*sqrt(N_1*N_2)*VDI(i))/M;
else
var_2=((Q_g_foundation*sqrt(N_1*N_2)*VDI(i))/M)+beta(1);
end

var_all=var_1+var_2^2;

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




