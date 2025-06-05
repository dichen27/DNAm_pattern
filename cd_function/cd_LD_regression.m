function [h_2_final,constant]=cd_LD_regression(R_2_moment,t_2_all,N_sample_size)



% lin


% m=2
% 
% R_2_moment=R2_all;
% t_2_all=t_all(:,m);

%% first regression of Nh^2/M
x=R_2_moment;
y=t_2_all;
% N_sample_size=size(cov_14,1)


x_one=[ones(size(x,1),1),x];


[b,bint,r,rint,stats] = regress(y,x_one);
% op_1=fc_ori_pat12-cov_pat12_one*b+b(1)*ones(size(fc_pat_zscore,1),1);%y-kx
% op_1=y-x_one*b;%r


beta_foundation=b(2);
%% h_foundation
M_num_voxel=length(R_2_moment);

h_2_foundation=(beta_foundation*M_num_voxel)/N_sample_size
%% Weighted
Weight=[];
for i=1:length(R_2_moment)
    
    var_lin=(1+beta_foundation*R_2_moment(i))^2;
    lin=1/var_lin;
    
    Weight=[Weight;lin];
end


%% beta_final
X=[ones(length(x),1),x];

% column to matrix
W=zeros(length(Weight),length(Weight));
for i=1:length(Weight)
W(i,i)=Weight(i);
end



lin=inv(X'*W*X);
beta_final=lin*X'*W*y;

% Res=y-X*beta_final;

beta_final=beta_final(2);

constant=beta_final(1);
%%
h_2_final=(beta_final*M_num_voxel)/N_sample_size;




end