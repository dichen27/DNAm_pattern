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