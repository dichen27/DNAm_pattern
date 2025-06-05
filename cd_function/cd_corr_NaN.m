function [r,p,n]=cd_corr_NaN(A,B)
% size of bb:  94 * 94

indx = ~(isnan(A) | isnan(B));
n=sum(indx);
[r,p]=corr(A(indx),B(indx));


end