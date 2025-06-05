function [a,b]=cd_bootstrap_95CI(r_bootstat)




a=mean(r_bootstat)-1.96*std(r_bootstat);
b=mean(r_bootstat)+1.96*std(r_bootstat);









end