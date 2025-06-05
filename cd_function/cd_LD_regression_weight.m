function [h_2_foundation,constant_foundation]=cd_LD_regression_weight(R_2_moment,t_attention,N_sample_size,weight_matrix)
% function  cd_LD_regression_weight(R_2_moment,t_attention,N_sample_size,weight_matrix)


% h_1=cd_LD_regression_weight(R2_all,t_all(:,m),N_sample_size,weight_R_matrix)

% lin
% R_2_moment=R2_all;
% t_2_all=t_all(:,m);
% N_sample_size
% weight_matrix=weight_R_matrix;

t_2_all=[];
for i=1:length(t_attention)
 t_2_all=[t_2_all;t_attention(i)^2];
end
%% first regression improve
M=length(R_2_moment);
x=R_2_moment;
y=t_2_all;

% W=zeros(M,M);
% for i=1:length(R_2_moment)
% W(i,i)=1/R_2_moment(i);% 1/L_j
% end

W=weight_matrix;

X=[ones(length(x),1),x];

lin=inv(X'*W*X);
beta=lin*X'*W*y;

% Res=y-X*beta_final;

beta_foundation=beta(2);

%% h_foundation
constant_foundation=beta(1);
h_2_foundation=(beta_foundation*M)/N_sample_size;

%% second regression
% 
% % %
% W=zeros(M,M);
% for i=1:length(R_2_moment)
%     
%     var_lin=(1+beta_foundation*R_2_moment(i))^2;
%     lin=1/var_lin; %1/var
% %     lin=1/sqrt(var_lin); %1/std
%         
%     W(i,i)=lin;
% end
% 
% 
% W=W*weight_matrix;
% 
% 
% % beta_final
% X=[ones(length(x),1),x];
% lin=inv(X'*W*X);
% beta=lin*X'*W*y;
% 
% % Res=y-X*beta;
% 
% beta_final=beta(2);
% constant=beta(1)
% 
% 
% %% Heritability
% 
% h_2_final=(beta_final*M)/N_sample_size


end