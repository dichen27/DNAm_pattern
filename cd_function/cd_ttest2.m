
function [t_improve, p_improve,df_improve]=cd_ttest2(brain_a_nor, grouplabel_a_nor, cov_a_nor)


cov_a_nor(:,find(sum(cov_a_nor)==0))=[];% dele the all 0 cov_a_nor in cov

t_improve=[];
p_improve=[];
df_improve=[];
for i=1:size(brain_a_nor,2)
    
%     i/size(brain_a_nor,2)
data_lin=[brain_a_nor(:,i),grouplabel_a_nor,cov_a_nor];

% dele NaN
data_lin(any(isnan(data_lin), 2),:) = [];

% ttest2_cov_improve
[tval_lin,p_lin] = ttest2_cov_improve(data_lin(:,1),data_lin(:,2),data_lin(:,3:end));

% df
df_lin=size(data_lin,1)-2-size(data_lin(:,3:end),2);


t_improve=[t_improve;tval_lin];
p_improve=[p_improve;p_lin];
df_improve=[df_improve;df_lin];
end

end

function [SSE, b] = regress_wei_improve(y,X)
% [b,r,SSE,SSR, T] = y_regress_ss(y,X)
% Perform regression.
% Revised from MATLAB's regress in order to speed up the calculation.
% Input:
%   y - Independent variable.
%   X - Dependent variable.
% Output:
%   b - beta of regression model.
%   r - residual.
%   SSE - The sum of squares of error.
%   SSR - The sum of squares of regression.
%   T - T value for each beta.


b = X \ y;
yhat = X*b;                     % Predicted responses at each data point.
r = y-yhat;                     % Residuals.
SSE=sum(r.^2);
end




function [tval, pval] = ttest2_cov_improve(DependentVariable, GroupLabel, Covariate)

index = sum(DependentVariable,2);
DependentVariable(isnan(index),:) = [];
GroupLabel(isnan(index)) = [];
Covariate(isnan(index),:) = [];

Df_E = size(DependentVariable,1) - 2 - size(Covariate,2);
SSE_H = regress_wei_improve(DependentVariable,[ones(size(DependentVariable,1),1),Covariate]);
SSE_H(isnan(SSE_H)) = 1;
% Calulate SSE
[SSE, b] = regress_wei_improve(DependentVariable,[ones(size(DependentVariable,1),1),GroupLabel,Covariate]);
SSE(isnan(SSE)) = 1;

% Calculate F
F = ((SSE_H-SSE)/1)./(SSE./Df_E);
pval =1-fcdf(F,1,Df_E);
tval = sqrt(F).*sign(b(2,:));
end