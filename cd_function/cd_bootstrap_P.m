function p_value=cd_bootstrap_P(Covariance,bootstat)






num_bootstrap=length(bootstat);


if Covariance>0
    p_value=length(find(bootstat<=0))/num_bootstrap;

else
p_value=length(find(bootstat>=0))/num_bootstrap;

end












end