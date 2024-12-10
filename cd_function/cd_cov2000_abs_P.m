function P=cd_cov2000_abs_P(cov)
% size of bb:  94 * 94

load('/share/inspurStorage/home1/chendi/Documents/LD_regression/MID_SST_mix/Permutation_final/Combine_old_2000.mat');

cov_abs=abs(cov);
P=length(find(cov_2000_abs>=cov_abs))/length(cov_2000_abs);

end